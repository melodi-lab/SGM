# SGM
Bipartite matching generalizations for peptide identification in tandem mass spectrometry


# User document SGM
  
Step 1: Run tide-index of <a href="http://crux.ms/">Crux</a> to generate decoy set.

CRUX tide-index --output-dir [OUTPUT] --clip-nterm-methionine T \
   --missed-cleavages [Missed cleavage] --peptide-list T \
   --mods-spec [MODS] --nterm-peptide-mods-spec [Nterm] \  
   --enzyme [ENZYME] \
   --overwrite T \ 
   [TARGETS] [tideindex] 
  

    See <a href="http://crux.ms/commands/tide-index.html">Crux</a> for desired parameters.

    Save tide-index.peptides.target.txt and tide-index.peptides.decoy.txt for input of next step. 

   

   Step 2: Run "write_msms_bipartite.py" to generate input file for SGM. 
 

  python $SHARDDATA \ 
    --spectra [MS2FILE] \ 
    --num_spectra $(grep -c '^S' [MS2FILE]) \
    --targets tide-index.peptides.target.txt \ 
    --decoys tide-index.peptides.decoy.txt \ 
    --tol [TOL] \ 
    --shards 1 \ 
    --output_dir [ENCODE] \ 
    --ppm \ 
    --charges [CHARGE] 

  for pickle in ENCODE/shard-*.pickle
  do 
    python write_msms_bipartite.py--output [ENCODE] --pickle $pickle \ 
    --max_mass 2000 --targets_per_spectrum 100000 \ 
    --normalize fastSequestTransform  \ 
    --bin_width 1.0005079 --bin_offset 0.68 --charge [CHARGE] \ 
    --bipartite_file [bipartite-OUTPUT] \
    --xcorr_ident_file [xcorrIdent-OUTPUT] \ 
    --do_not_write_all_psms --foreground
  done

  [ENCODE] is the output folder. [TOL] is the MS1 torrlerance in ppm. [CHARGE] is the precursor charge. [bipartite-OUTPUT] and [xcorrIdent-OUTPUT] are two output files.

 
Step 3: Run SGM. 
 Complie and run with following parameters.
 
option: -h  show help information
      -c  charge <br />
     -i  [bipartite-OUTPUT] 
      -o  output file
       -m  ms2file 

          </P>

Optional: Merge different charge state.
run merge.py

Optional: Draw absulte plot.
run absplot.py


</div>
  


