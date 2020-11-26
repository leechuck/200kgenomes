new java.util.zip.GZIPInputStream(new FileInputStream("/encrypted/e3001/newexomes/200k-chr1.vcf.gz")).splitEachLine("\t") { line ->
    if (line[0].startsWith("#CHROM")) {
	line[9..-1].each {
	    println it
	}
	System.exit(0)
    }
}
