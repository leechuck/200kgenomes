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

// TODO: some stuff to get gene and individual number
// count samples and genes
//def samples = new LinkedHashSet()
def genes = new LinkedHashSet()
l.each { chr ->
    new java.util.zip.GZIPInputStream(new FileInputStream("/ibex/scratch/projects/c2014/robert/chr${chr}-quantiles.txt.gz")).splitEachLine("\t") { line ->
	def gene = line[0]
//	def sample = line[1]
	genes.add(gene)
//	samples.add(samples)
    }
}
PrintWriter fout = new PrintWriter(new BufferedWriter(new FileWriter("/ibex/scratch/projects/c2014/robert/genes.txt")))
genes.eachWithIndex { gene, idx -> fout.println("$gene\t$idx") }
fout.flush()
fout.close()
fout = new PrintWriter(new BufferedWriter(new FileWriter("/ibex/scratch/projects/c2014/robert/samples.txt")))
samples.eachWithIndex { sample, idx -> fout.println("$sample\t$idx") }
fout.flush()
fout.close()
System.exit(1)

DoubleFactory2D factory = DoubleFactory2D.dense
DoubleMatrix2D hetmat = factory.make(200000, 20000) // this is an individual x gene matrix, about 200,000 x 20,000 double values -> 32 GB
DoubleMatrix2D hommat = factory.make(200000, 20000) // this is an individual x gene matrix, about 200,000 x 20,000 double values -> 32 GB

//def cmap = [:].withDefault { [:].withDefault { [:] } } // chr -> pos -> alt -> val
l.each { chr ->
    new java.util.zip.GZIPInputStream(new FileInputStream("/ibex/scratch/projects/c2014/robert/chr${chr}-quantiles.txt.gz")).splitEachLine("\t") { line ->    
	def gene = line[0]
	def pid = line[1]
	def hetval = new Float(line[2])
	def hetrank = new Float(line[6]) / new Float(line[7]) // kind of the quantile, mean rank / list size; interpolated rank
	def homval = new Float(line[8])
	def homrank = new Float(line[12]) / new Float(line[13])
    }
}
