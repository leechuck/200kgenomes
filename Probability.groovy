// reads the values, generates matrix of patient x gene, then computes the Friedman statistic; then performs Tukey's honest significance test

@Grab(group='org.codehaus.gpars', module='gpars', version='1.0.0')
  @Grab(group='colt', module='colt', version='1.2.0')

import static groovyx.gpars.GParsPool.withPool
import cern.colt.list.*
import cern.colt.matrix.*
import cern.jet.stat.*

def l = []
(1..22).each { l << it }
l << "X"
l << "Y"

def samples = [:]
def count = 0
new File("/ibex/scratch/projects/c2014/robert/samples.txt").eachLine { line ->
    samples[line.trim()] = count
    count += 1
}
count = 0
def genes = [:]
new File("/ibex/scratch/projects/c2014/robert/genes.txt").eachLine { line ->
    genes[line.trim()] = count
    count += 1
}

def filenames = []
def process = ["sh", "-c", "ls /ibex/scratch/projects/c2014/robert/chr*-*-quantiles.txt.gz"].execute()
process.waitFor()
process.text.eachLine { fn ->
    filenames << fn
}

Float[][] hetmat = new Float[samples.size()][genes.size()] // this is a sample x gene matrix, about 200,000 x 20,000 double values -> 32 GB
Float[][] hommat = new Float[samples.size()][genes.size()] // this is an individual x gene matrix, about 200,000 x 20,000 double values -> 32 GB

//def cmap = [:].withDefault { [:].withDefault { [:] } } // chr -> pos -> alt -> val
filenames.each { fn ->
    println fn
    new java.util.zip.GZIPInputStream(new FileInputStream(fn)).splitEachLine("\t") { line ->    
	def gene = line[0]
	def pid = line[1]
	def hetval = new Float(line[2])
	def hetrank = new Float(line[6]) / new Float(line[7]) // kind of the quantile, mean rank / list size; interpolated rank
	def homval = new Float(line[8])
	def homrank = new Float(line[12]) / new Float(line[13])
	def gpos = genes[gene]
	def ppos = samples[pid]
	hetmat[ppos][gpos] = hetval
	hommat[ppos][gpos] = homval
    }
}
