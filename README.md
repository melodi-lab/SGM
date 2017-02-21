

<HTML>
   <HEAD>
      <TITLE>
        
      </TITLE>
    <!-- Bootstrap Core CSS -->
	<link rel="stylesheet" href="assets/css/bootstrap.css" rel="stylesheet">
	<!-- Template CSS -->
	<link rel="stylesheet" href="assets/css/animate.css" rel="stylesheet">
	<link rel="stylesheet" href="assets/css/font-awesome.css" rel="stylesheet">
	<link rel="stylesheet" href="assets/css/nexus.css" rel="stylesheet">
	<link rel="stylesheet" href="assets/css/responsive.css" rel="stylesheet">
	<link rel="stylesheet" href="assets/css/custom.css" rel="stylesheet">
	<!-- Google Fonts-->
	<link href="http://fonts.googleapis.com/css?family=Lato:400,300" rel="stylesheet" type="text/css">
	<link href="http://fonts.googleapis.com/css?family=Open+Sans:400,300" rel="stylesheet" type="text/css">
    
   </HEAD>
<BODY>
<div id="content" class="container">

<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=AM_CHTML"></script>


  <H1> Bipartite matching generalizations for peptide identification in tandem mass spectrometry</H1>
  <H2> How to use SGM </H2>
  
     <H3>Step 1: Run tide-index of <a href="http://crux.ms/">Crux</a> to generate decoy set.</H3>
     <P>
    CRUX tide-index --output-dir [OUTPUT] --clip-nterm-methionine T \<br /> 
    --missed-cleavages [Missed cleavage] --peptide-list T \ <br />
    --mods-spec [MODS] --nterm-peptide-mods-spec [Nterm] \ <br />    
    --enzyme [ENZYME] \ <br />
    --overwrite T \ <br />
    [TARGETS] [tideindex] <br />
    <br />

    See <a href="http://crux.ms/commands/tide-index.html">Crux</a> for desired parameters. <br />

    Save tide-index.peptides.target.txt and tide-index.peptides.decoy.txt for input of next step. <br />

  </P>
   

   <H3>Step 2: Run "write_msms_bipartite.py" to generate input file for SGM. </H3>
  <P>

  python $SHARDDATA \ <br />
    --spectra [MS2FILE] \ <br /> 
    --num_spectra $(grep -c '^S' [MS2FILE]) \ <br />
    --targets tide-index.peptides.target.txt \ <br />
    --decoys tide-index.peptides.decoy.txt \ <br />
    --tol [TOL] \ <br />
    --shards 1 \ <br />
    --output_dir [ENCODE] \ <br />
    --ppm \ <br />
    --charges [CHARGE] <br /> <br />

  for pickle in ENCODE/shard-*.pickle <br />
  do <br />
    python write_msms_bipartite.py--output [ENCODE] --pickle $pickle \ <br />
    --max_mass 2000 --targets_per_spectrum 100000 \ <br />
    --normalize fastSequestTransform  \ <br />
    --bin_width 1.0005079 --bin_offset 0.68 --charge [CHARGE] \ <br />
    --bipartite_file [bipartite-OUTPUT] \ <br />
    --xcorr_ident_file [xcorrIdent-OUTPUT] \ <br />
    --do_not_write_all_psms --foreground <br />
  done <br />

  [ENCODE] is the output folder. [TOL] is the MS1 torrlerance in ppm. [CHARGE] is the precursor charge. [bipartite-OUTPUT] and [xcorrIdent-OUTPUT] are two output files.

  </P>
 <H3>Step 3: Run SGM. </H3>
 Complie and run with following parameters.
 </P>
option: -h  show help information <br />
<span style="display:inline-block; width: 55;"></span>        -c  charge <br />
<span style="display:inline-block; width: 55;"></span>          -i  [bipartite-OUTPUT] <br />
<span style="display:inline-block; width: 55;"></span>          -o  output file <br />
<span style="display:inline-block; width: 55;"></span>          -m  ms2file <br />

          </P>

<H3>Optional: Merge different charge state. </H3>
<P>run merge.py</P>

<H4>Optional: Draw absulte plot. </H3>
<P>run absplot.py</P>
<H2> TODO </H2>

<P>Make the codes more user friendly. </P>


</div>
</BODY>
</HTML>


