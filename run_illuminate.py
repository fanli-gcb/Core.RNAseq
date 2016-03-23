from illuminate import InteropDataset
import argparse
import sys

parser = argparse.ArgumentParser(description='Parse Illumina InterOp metrics.')
parser.add_argument('run', metavar='RUN_FOLDER', type=str, help='Path to Illumina Run folder')
parser.add_argument('--output', dest='output', help='Output filename (default is STDOUT)')
args = parser.parse_args()
myDataset = InteropDataset(args.run)
doError = True

if args.output:
	fout = open(args.output, 'w')
else:
	fout = sys.stdout

# metadata
metadata = myDataset.Metadata()
fout.write("Experiment_Name\t%s\n" % metadata.experiment_name)
fout.write("Start_Datetime\t%s\n" % metadata.start_datetime)
fout.write("Chemistry\t%s\n" % metadata.chemistry)

# tile metrics
tilemetrics = myDataset.TileMetrics()
fout.write("Mean_Cluster_Density\t%i\n" % tilemetrics.mean_cluster_density)
fout.write("Mean_PF_Cluster_Density\t%i\n" % tilemetrics.mean_cluster_density_pf)
fout.write("Total_Clusters\t%i\n" % tilemetrics.num_clusters)
fout.write("Total_PF_Clusters\t%i\n" % tilemetrics.num_clusters_pf)
fout.write("Percent_PF_Clusters\t%f\n" % tilemetrics.percent_pf_clusters)
for read_num in range(tilemetrics.num_reads):
	fout.write("Read%i_Phasing\t%f\n" % (read_num+1, tilemetrics.mean_phasing[read_num]*100))
	fout.write("Read%i_Prephasing\t%f\n" % (read_num+1, tilemetrics.mean_prephasing[read_num]*100))
	fout.write("Read%i_Aligned\t%f\n" % (read_num+1, tilemetrics.mean_aligned[read_num]))
# quality metrics
qualitymetrics = myDataset.QualityMetrics()
#for key in dir(qualitymetrics):
#	fout.write("%s\n" % key)
for read in qualitymetrics.read_config:
	fout.write("Read%i_Percent_Q30\t%f\n" % (read['read_num'], qualitymetrics.get_qscore_percentage(30, read['read_num']-1)))
fout.write("Overall_Percent_Q30\t%f\n" % qualitymetrics.get_qscore_percentage(30, -1))
# error metrics
try:
	errormetrics = myDataset.ErrorMetrics()
except:
	doError = False
if doError:
#	df = errormetrics.df
#	total_cycles = 0
#	df['read'] = 1
#	for read in errormetrics.read_config:
#		for i in range(1, read['cycles']):
#			df.loc[df['cycle']==(i+total_cycles), 'read'] = read['read_num']
#		total_cycles += read['cycles']
#	errormetrics.mean_error_rate = df.pivot_table('rate', 'read', aggfunc="mean")
	for read in errormetrics.read_config:
		if not read['is_index']:
			fout.write("Read%i_Error_Rate\t%f\n" % (read['read_num'], errormetrics.mean_error_rate[read['read_num']]))

fout.close()

##qualitymetrics = myDataset.QualityMetrics()
##indexmetrics = myDataset.IndexMetrics()
##controlmetrics = myDataset.ControlMetrics()
##corintmetrics = myDataset.CorrectedIntensityMetrics()
##extractionmetrics = myDataset.ExtractionMetrics()
##errormetrics = myDataset.ErrorMetrics()


