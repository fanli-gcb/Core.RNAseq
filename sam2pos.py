import sys, optparse, itertools, warnings, traceback, os.path

import HTSeq

from collections import defaultdict

class UnknownChrom( Exception ):
   pass
   
def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

def count_reads_in_features( sam_filename, gff_filename, stranded, 
      overlap_mode, feature_type, id_attribute, quiet, minaqual, samout, allow_ambiguous, allow_nonunique ):
      
   def write_to_samout( r, assignment ):
      if samoutfile is None:
         return
      if not pe_mode:
         r = (r,)
      for read in r:
         if read is not None:
            samoutfile.write( read.original_sam_line.rstrip() + 
               "\tXF:Z:" + assignment + "\n" )
      
   if quiet:
      warnings.filterwarnings( action="ignore", module="HTSeq" ) 
      
   if samout != "":
      samoutfile = open( samout, "w" )
   else:
      samoutfile = None
      
   features = HTSeq.GenomicArrayOfSets( "auto", stranded != "no" )
   features_dict = defaultdict(list)
   counts = {}

   # Try to open samfile to fail early in case it is not there
   if sam_filename != "-":
      open( sam_filename ).close()
      
   gff = HTSeq.GFF_Reader( gff_filename )   
   i = 0
   try:
      for f in gff:
         if f.type == feature_type:
            try:
               feature_id = f.attr[ id_attribute ]
            except KeyError:
               sys.exit( "Feature %s does not contain a '%s' attribute" % 
                  ( f.name, id_attribute ) )
            if stranded != "no" and f.iv.strand == ".":
               sys.exit( "Feature %s at %s does not have strand information but you are "
                  "running htseq-count in stranded mode. Use '--stranded=no'." % 
                  ( f.name, f.iv ) )
            features[ f.iv ] += feature_id
            counts[ f.attr[ id_attribute ] ] = 0
            features_dict[ f.attr[ id_attribute ] ].append(f)
         i += 1
         if i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d GFF lines processed.\n" % i )
   except:
      sys.stderr.write( "Error occured in %s.\n" % gff.get_line_number_string() )
      raise
      
   if not quiet:
      sys.stderr.write( "%d GFF lines processed.\n" % i )
   
   sys.stderr.write( "Sorting exons from GFF file.\n" )
   for key, value in features_dict.items():
   		if features_dict[key][0].iv.strand == "-":
   				features_dict[key] = sorted(features_dict[key], key=lambda feat: feat.iv.start, reverse=True)
   		else:
   				features_dict[key] = sorted(features_dict[key], key=lambda feat: feat.iv.start, reverse=False)
      
   if len( counts ) == 0 and not quiet:
      sys.stderr.write( "Warning: No features of type '%s' found.\n" % feature_type )
   
   try:
      if sam_filename != "-":
         read_seq = HTSeq.SAM_Reader( sam_filename )
         first_read = iter(read_seq).next()
      else:
         read_seq = iter( HTSeq.SAM_Reader( sys.stdin ) )
         first_read = read_seq.next()
         read_seq = itertools.chain( [ first_read ], read_seq )
      pe_mode = first_read.paired_end
   except:
      sys.stderr.write( "Error occured when reading first line of sam file.\n" )
      raise

   try:
      if pe_mode:
         read_seq_pe_file = read_seq
         read_seq = HTSeq.pair_SAM_alignments( read_seq )
      empty = 0
      ambiguous = 0
      notaligned = 0
      lowqual = 0
      nonunique = 0
      i = 0   
      for r in read_seq:
         i += 1
         if not pe_mode:
            if not r.aligned:
               notaligned += 1
               write_to_samout( r, "not_aligned" )
               continue
            try:
               if ((allow_nonunique == "no") and (r.optional_field( "NH" ) > 1)):
                  write_to_samout( r, "alignment_not_unique" )
                  nonunique += 1
                  continue
            except KeyError:
               pass
            if r.aQual < minaqual:
               lowqual += 1
               write_to_samout( r, "too_low_aQual" )
               continue
            if stranded != "reverse":
               iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" )
            else:
               iv_seq = ( invert_strand( co.ref_iv ) for co in r.cigar if co.type == "M" )            
         else:
            if r[0] is not None and r[0].aligned:
               if stranded != "reverse":
                  iv_seq = ( co.ref_iv for co in r[0].cigar if co.type == "M" )
               else:
                  iv_seq = ( invert_strand( co.ref_iv ) for co in r[0].cigar if co.type == "M" )
            else:
               iv_seq = tuple()
            if r[1] is not None and r[1].aligned:            
               if stranded != "reverse":
                  iv_seq = itertools.chain( iv_seq, 
                     ( invert_strand( co.ref_iv ) for co in r[1].cigar if co.type == "M" ) )
               else:
                  iv_seq = itertools.chain( iv_seq, 
                     ( co.ref_iv for co in r[1].cigar if co.type == "M" ) )
            else:
               if ( r[0] is None ) or not ( r[0].aligned ):
                  write_to_samout( r, "not_aligned" )
                  notaligned += 1
                  continue         
            try:
               if (allow_nonunique == "no") and (( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 )):
                  nonunique += 1
                  write_to_samout( r, "alignment_not_unique" )
                  continue
            except KeyError:
               pass
            if ( r[0] and r[0].aQual < minaqual ) or ( r[1] and r[1].aQual < minaqual ):
               lowqual += 1
               write_to_samout( r, "too_low_aQual" )
               continue         
         
         try:
            iv_seq = list(iv_seq)
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     fs = fs.union( fs2 )
            elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
               fs = None
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     if len(fs2) > 0 or overlap_mode == "intersection-strict":
                        if fs is None:
                           fs = fs2.copy()
                        else:
                           fs = fs.intersection( fs2 )
            else:
               sys.exit( "Illegal overlap mode." )
#            sys.stderr.write( "fs = %s with len = %d allow_ambiguous=%s\n" % (fs, len(fs), allow_ambiguous) )
            if fs is None or len( fs ) == 0:
               write_to_samout( r, "no_feature" )
               empty += 1
            elif ((len( fs ) > 1) and (allow_ambiguous == "no")):
               write_to_samout( r, "ambiguous[" + '+'.join( fs ) + "]" )
               ambiguous += 1
            else:
#               write_to_samout( r, list(fs)[0] )
#               sys.stderr.write( "iv_seq is now %s\n" % iv_seq )
#               sys.exit ("stopping for now" )
               for iv in iv_seq:
#               		sys.stderr.write( "iv in iv_seq is %s and getting %s\n" % (iv, features[iv]) )
               		for iv2, fs2 in features[ iv ].steps():
#						         sys.stderr.write( "iv2 = %s and fs2 = %s\n" % (iv2, fs2) )
						         if len(fs2) == 0:
						         		continue
						         for fsi in fs2:
#						         		sys.stderr.write( "fsi = %s\n" % fsi )
						         		offset = 0
#						         		write_to_samout( r, fsi[0] )
						         		for exon in features_dict[ fsi ]:
						         				if ((iv2.start >= exon.iv.start) and (iv2.end <= exon.iv.end)):
#						         						sys.stderr.write("found matching exon %s\n" % exon)
						         						if (exon.iv.strand == "+"):
						         								offset += (iv2.start - exon.iv.start)
						         						else:
						         								offset += (exon.iv.end - iv2.end)
#						         								sys.stderr.write("matching exon new offset %d\n" % offset)
						         						break
						         				else:
						         						offset += (exon.iv.end - exon.iv.start)
#						         						sys.stderr.write( "skipping exon %s - new offset %d\n" % (exon, offset) )
#						         				sys.stderr.write( "mapping read %s to offset %d\n" % (r, offset) )
#						         				sys.stderr.write( "dir(r) is %s" % "\t".join(list(dir(r))) )
						         		if pe_mode:
						         			rname = r[0].read.name if r[0] is not None else r[1].read.name
						         		else:
						         			rname = r.read.name
						         		print "%s\t%d\t%d\t%s" % (fsi, offset, offset + (iv2.end-iv2.start-1), rname) # output is 0-based, inclusive on both ends
         except UnknownChrom:
            if not pe_mode:
               rr = r 
            else: 
               rr = r[0] if r[0] is not None else r[1]
            empty += 1
            #if not quiet:
            #   sys.stderr.write( ( "Warning: Skipping read '%s', because chromosome " +
            #      "'%s', to which it has been aligned, did not appear in the GFF file.\n" ) % 
            #      ( rr.read.name, iv.chrom ) )

         if i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d sam %s processed.\n" % ( i, "lines " if not pe_mode else "line pairs" ) )

   except:
      if not pe_mode:
         sys.stderr.write( "Error occured in %s.\n" % read_seq.get_line_number_string() )
      else:
         sys.stderr.write( "Error occured in %s.\n" % read_seq_pe_file.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d sam %s processed.\n" % ( i, "lines " if not pe_mode else "line pairs" ) )
         
   if samoutfile is not None:
      samoutfile.close()

#   for fn in sorted( counts.keys() ):
#      print "%s\t%d" % ( fn, counts[fn] )
#   print "no_feature\t%d" % empty
#   print "ambiguous\t%d" % ambiguous
#   print "too_low_aQual\t%d" % lowqual
#   print "not_aligned\t%d" % notaligned
#   print "alignment_not_unique\t%d" % nonunique

      
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] sam_file gff_file",
      
      description=
         "This script takes an alignment file in SAM format and a " +
         "feature file in GFF format and calculates for each feature " +
         "the number of reads mapping to it. See " +
         "http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for details.",
         
      epilog = 
         "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
         "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
         "Public License v3. Part of the 'HTSeq' framework, version %s." % HTSeq.__version__ )
         
   optParser.add_option( "-m", "--mode", type="choice", dest="mode",
      choices = ( "union", "intersection-strict", "intersection-nonempty" ), 
      default = "union", help = "mode to handle reads overlapping more than one feature" +
         "(choices: union, intersection-strict, intersection-nonempty; default: union)" )
         
   optParser.add_option( "-s", "--stranded", type="choice", dest="stranded",
      choices = ( "yes", "no", "reverse" ), default = "yes",
      help = "whether the data is from a strand-specific assay. Specify 'yes', " +
         "'no', or 'reverse' (default: yes). " +
         "'reverse' means 'yes' with reversed strand interpretation" )    
         
   optParser.add_option( "-g", "--allow_ambiguous", type="choice", dest="allow_ambiguous",
   		choices = ( "yes", "no" ), default = "no",
   		help = "whether to allow mappings to multiple transcripts" )
   		
   optParser.add_option( "-n", "--allow_nonunique", type="choice", dest="allow_nonunique",
   		choices = ( "yes", "no" ), default = "no",
   		help = "whether to allow nonunique mapping" )
      
   optParser.add_option( "-a", "--minaqual", type="int", dest="minaqual",
      default = 0,
      help = "skip all reads with alignment quality lower than the given " +
         "minimum value (default: 0)" )
      
   optParser.add_option( "-t", "--type", type="string", dest="featuretype",
      default = "exon", help = "feature type (3rd column in GFF file) to be used, " +
         "all features of other type are ignored (default, suitable for Ensembl " +
         "GTF files: exon)" )
         
   optParser.add_option( "-i", "--idattr", type="string", dest="idattr",
      default = "gene_id", help = "GFF attribute to be used as feature ID (default, " +
      "suitable for Ensembl GTF files: gene_id)" )

   optParser.add_option( "-o", "--samout", type="string", dest="samout",
      default = "", help = "write out all SAM alignment records into an output " +
      "SAM file called SAMOUT, annotating each line with its feature assignment " +
      "(as an optional field with tag 'XF')" )

   optParser.add_option( "-q", "--quiet", action="store_true", dest="quiet",
      help = "suppress progress report and warnings" )

   if len( sys.argv ) == 1:
      optParser.print_help()
      sys.exit(1)

   (opts, args) = optParser.parse_args()

   if len( args ) != 2:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide two arguments.\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
   warnings.showwarning = my_showwarning
   try:
      count_reads_in_features( args[0], args[1], opts.stranded, 
         opts.mode, opts.featuretype, opts.idattr, opts.quiet, opts.minaqual,
         opts.samout, opts.allow_ambiguous, opts.allow_nonunique )
   except:
      sys.stderr.write( "Error: %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )

def my_showwarning( message, category, filename, lineno = None, line = None ):
   sys.stderr.write( "Warning: %s\n" % message )

if __name__ == "__main__":
   main()

