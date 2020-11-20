def varset = [:] // .withDefault{ [:].withDefault { [:] } } // chr -> pos -> alt -> score
def col2patient = [:]
def count = 0

def chr = args[0].substring(lastIndexOf("/")+1).replaceAll("200k-chr","").replaceAll(".vcf.gz","")
new java.util.zip.GZIPInputStream(new FileInputStream(args[0])).splitEachLine("\t") { line ->
    if (line[0].startsWith("#CHROM")) {
	line[9..-1].eachWithIndex { id, idx ->
	    col2patient[idx] = id
	}
    } else if (! line[0].startsWith("#")) {
	
	def chr = new Integer(line[0])
	def pos = new Integer(line[1])
	def alt = line[4]
	if (varset[chr] == null) {
	    varset[chr] = [:]
	}
	if (varset[chr][pos] == null) {
	    varset[chr][pos] = [:]
	}
	varset[chr][pos][alt] = new Float(0.0)
    }
    count += 1
    if (count % 1000 == 0) {
	println count + " variants processed"
    }
}

PrintWriter fout = new PrintWriter(new BufferedWriter(new FileWriter("/ibex/scratch/projects/c2014/robert/cadd-scores-${chr}.txt")))

count = 0
new java.util.zip.GZIPInputStream(new FileInputStream("/ibex/sw/csi/cadd/1.6/el7.7_miniconda2/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/whole_genome_SNVs.tsv.gz")).splitEachLine("\t") { line ->
    if (! line[0].startsWith("#")) {
	def chr = new Integer(line[0])
	def pos = new Integer(line[1])
	def alt = line[3]
	if (varset[chr] && varset[chr][pos] && varset[chr][pos][alt]) {
	    // varset[chr][pos][alt] = new Float(line[5])
	    fout.println("$chr\t$pos\t$alt\t" + line[5])
	}
    }
    count += 1
    if (count % 1000000 == 1) {
	println count + " CADD variants processed"
    }
}
fout.flush()
fout.close()
