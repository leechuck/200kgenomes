
// reads the CADD annotations; fills in gene values for regulatory regions

def coordinates = [:]
new java.util.zip.GZIPInputStream(new FileInputStream("/ibex/scratch/projects/c2014/robert/Homo_sapiens.GRCh38.101.chr.gtf.gz")).splitEachLine("\t") { line ->
    if (! line[0].startsWith("#")) {
	if (line[2] == "gene") {
	    def chr = line[0]
	    if (coordinates[chr] == null) {
		coordinates[chr] = []
	    }
	    def start = new Integer(line[3])
	    def end = new Integer(line[4])
	    def info = line[8].split(";")
	    def geneid = info[0].replaceAll("gene_id \"","").replaceAll("\"","")
	    Expando exp = new Expando()
	    exp.chr = chr
	    exp.start = start
	    exp.end = end
	    exp.geneid = geneid
	    coordinates[chr] << exp
	}
    }
}

coordinates.each { chr, l ->
    for (def i = 0 ; i < l.size() ; i++) {
	Expando exp = l[i]
	if (i == 0) {
	    exp.start = exp.start - 10000 // we start 1000 bp before in case of the first gene
	} else if (i == l.size() - 1) { // we end 1000 bp after in case of the last gene
	    exp.end = exp.end + 10000
	} else {
	    Expando exp2 = l[i-1] // get the gene before
	    def middle = Math.ceil((exp.start + exp2.end)/2).toInteger()
	    if (exp2.end < middle-1) {
		exp2.end = middle - 1
	    }
	    if (exp.start > middle) {
		exp.start = middle
	    }
	}
    }
}

def count = 0
new java.util.zip.GZIPInputStream(new FileInputStream(args[0])).splitEachLine("\t") { line ->
    if (! line[0].startsWith("#")) {
	def chr = line[0]
	def pos = new Integer(line[1])
	def gene = coordinates[chr][count]
	while (pos > gene.end) {
	    count += 1
	    gene = coordinates[chr][count]
	}
	def ref = line[2]
	def alt = line[3]
	def score = new Float(line[-1])
	println "$chr\t$pos\t$ref\t$alt\t$score\t${gene.geneid}"
    }
}
