@Grab(group='org.codehaus.gpars', module='gpars', version='1.0.0')
  @Grab(group='colt', module='colt', version='1.2.0')

// computed probability of having a variant higher than x in cohort
import static groovyx.gpars.GParsPool.withPool
import cern.colt.list.*
import cern.jet.stat.*

// read CADD values in map
def cadd = [:]
def pos2gene = [:]
new java.util.zip.GZIPInputStream(new FileInputStream(args[0])).splitEachLine("\t") { line ->
    def chr = line[0]
    if (cadd[chr] == null) {
	cadd[chr] = [:]
    }
    def pos = new Integer(line[1])
    if (cadd[chr][pos] == null) {
	cadd[chr][pos] = [:]
    }
    def ref = line[2]
    def alt = line[3]
    def score = new Float(line[4])
    cadd[chr][pos][alt] = score
    def gene = line[5]
    pos2gene[pos] = gene
}

def hets = [:].withDefault { 0 } // keeps max for current gene only; index -> max val
def homs = [:].withDefault { 0 } // keeps max for current gene only; index -> max val
def current = ""
//PrintWriter fout = new PrintWriter(new BufferedWriter(new FileWriter(args[2])))
EPSILON = 0.001
def start = true
def index2id = [:]
withPool(16) {
    new java.util.zip.GZIPInputStream(new FileInputStream(args[1])).splitEachLine("\t") { line ->
	if (start && line[0].startsWith("#CHROM")) {
	    start = false
	    line[9..-1].eachWithIndex { id, idx ->
		index2id[idx] = id
	    }
	}
	if (! line[0].startsWith("#")) {
	    def chr = line[0]
	    def pos = new Integer(line[1])
	    def alt = line[4]
	    def score = cadd[chr][pos][alt]
	    def gene = pos2gene[pos]
	    if (gene != current) {

		DoubleArrayList l = new DoubleArrayList()
		DoubleArrayList l2 = new DoubleArrayList()
		hets.each { k, v -> l.add(v) }
		homs.each { k, v -> l2.add(v) }
		if (l.size()>0) {
		    l.sort()
		    l2.sort()
		    hets.each { k, v -> // hets include the homs, so we just iterate through the hets
			def hh = homs[k]
			def q = Descriptive.quantileInverse(l, v - EPSILON) // strictly LESS
			def hq = Descriptive.quantileInverse(l2, hh - EPSILON) // strictly LESS
			println "$gene\t${index2id[k]}\t$v\t$q\t$hh\t$hq"
		    }
		}
		// reset to new gene
		hets = [:].withDefault { 0 }
		homs = [:].withDefault { 0 }
		current = gene
	    }
	    line[9..-1].eachWithIndexParallel { geno, idx ->
		if (geno == "0/1") {
		    if (hets[idx] <= score) {
			hets[idx] = score
		    }
		}
		if (geno == "1/1") {
		    if (homs[idx] <= score) {
			homs[idx] = score
			}
			// add homs to the hets
		    if (hets[idx] <= score) {
			hets[idx] = score
		    }
		}
	    }
	}
    }
}

