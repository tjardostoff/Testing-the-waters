###############################################################################
# script name   : asv_table_long.awk
# description   : Return an ASV/OTU table in long format whilst subsetting it.
# usage         : to extract and rename samples listed in "corresp.tsv" from the 
#                 ASV table "my_table.tsv" and to only keep ASVs with at least
#                 10 reads and fount in at least 2 sample (amongst the sample 
#                 selection), run in a terminal:
#
#                 awk \
#                   -v mintotal=10 \
#                   -v minspread=2 \
#                   -v asv_id_index=1 \
#                   -f asv_table_longer.awk \
#                   corresp.tsv \
#                   my_table.tsv
#
###############################################################################

BEGIN {

  FS = "\t"
  OFS = "\t"
  print "asv_id","sample","nreads"

}

# list samples to keep from corresp.tsv

NR == FNR {

  samples_to_select[$1] = $2
  next

}

# list columns to keep in my_table.tsv and their new names

FNR == 1 {

  for (i = 1; i <= NF; ++i){

    if(samples_to_select[$i]) col_to_keep[i] = samples_to_select[$i]

  }

}

FNR > 1 {

  total = 0
  # to empty final_sample_summed_values
  split("", summed_values)

  for(i in col_to_keep){
    total += $i
  } 

  if(total >= mintotal){
    
    # sum value for each new sample name
    for(i in col_to_keep){

        if($i>0) {
          summed_values[col_to_keep[i]] += $i
        }

    }

    # only print ASV dtected in at least minspread samples 
    if(length(summed_values) >= minspread) {

      for(new_sample in summed_values){

        print $asv_id_index, new_sample, summed_values[new_sample]

      }

    }

  }

}