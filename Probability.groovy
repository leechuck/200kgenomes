// reads the values, generates matrix of patient x gene, then computes the Friedman statistic

@Grab(group='org.codehaus.gpars', module='gpars', version='1.0.0')
  @Grab(group='colt', module='colt', version='1.2.0')

import static groovyx.gpars.GParsPool.withPool
import cern.colt.list.*
import cern.colt.matrix.*
import cern.jet.stat.*

DoubleFactory2D factory = DoubleFactory2D.dense

DoubleMatrix2D mat = factory.make(5,5)

Double friedman(DoubleMatrix2D matrix) {
    
}

def l = []
(1..22).each { l << it }
l << "X"
l << "Y"

// read the CADD scores
def cmap = [:].withDefault { [:].withDefault { [:] } } // chr -> pos -> alt -> val
l.each { chr ->
    new java.util.zip.GZIPInputStream(new FileInputStream("/encrypted/e3001/newexomes/200k-chr${chr}.cadd.gz")).splitEachLine("\t") { line ->    
	def chr = line[0]
	def pos = new Integer(line[1])
	def alt = line[3]
	def score = new Float(line[4])
	def gene = line[5]
    }
}
