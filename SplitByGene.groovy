// splits processed files so they contain exactly XXX genes

def current = ""
def flag = true
def count = 0
def fcount = 0 // file counter
def XXX = 100
PrintWriter fout = null
new java.util.zip.GZIPInputStream(new FileInputStream(args[0])).eachLine { line ->
    if (flag) { // first
	def toks = line.split("\t")
	current = toks[5]
	flag = false
	fout = new PrintWriter(new BufferedWriter(new FileWriter(args[0]+"."+fcount)))
    } else {
	def toks = line.split("\t")
	def gene = toks[5]
	if (gene != current) {
	    count += 1
	    current = gene
	    if (count % XXX == 0) {
		fout.flush()
		fout.close()
		("gzip "+args[0]+"."+fcount).execute()
		fcount += 1
		fout = new PrintWriter(new BufferedWriter(new FileWriter(args[0]+"."+fcount)))
	    }
	    fout.println(line)
	} else {
	    fout.println(line)
	}
    }
    
}
("gzip "+args[0]+"."+fcount).execute()
fout.flush()
fout.close()
