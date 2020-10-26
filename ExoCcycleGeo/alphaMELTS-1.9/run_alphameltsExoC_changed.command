#!/usr/bin/perl -w

# run_alphamelts.command 1.9
use strict;

my (@argv2, $in_file, $melts_file, $out_file, $log_file, $batch_file, $column_file, $table_file, $table_row);
my ($path, $run_path, $batch, $subsol, %com, $oxide_list, @phases, @used_oxides, %all_oxides, @all_pt, @all_masses);
my ($windows, $program, $rs, $run, $delim, $line, @lines, @oldlines, $title);

$in_file = '/Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/ExoC/Lith_bndry_env.txt';

@argv2 = ();
until (@argv2) {

    if (@ARGV) {

 	@argv2 = @ARGV;

    }
    else {

	my @temp;

	print ("No command line switches specified! ");
	print ("Please enter switches now\n(or press return for default '-f $in_file').\n");

	$_ = <STDIN>;
	chomp;

	s/\s+$//; # Trim any trailing white space
	s/\'/\"/g; # Windows can only deal with name enclosed in double quotes; *nix can take single or double

	s/\s+-/ --/g; # Assume that filenames etc. cannot begin with '-'
	@temp = split /\s+-/;
	map { push @argv2, (split /\s+/, $_, 2); } @temp;

    }

    # -h prints the following help before quitting
    if (@argv2 && $argv2[0] eq '-h') {
	warn ("\n usage: run_alphamelts\.command \n\n".
	      " [-h] print this brief help\n".
	      " [-f settings_file] input file for environment variables etc. \n".
	      "              (default $in_file) \n".
	      " [-b batch_file] run in batch mode using the file given as command line input  \n".
	      " [-a] run in automatic mode (writes auto_batch.txt and uses it in batch mode) \n".
	      " [-m melts_file] file for initial compositions etc. \n".
	      " [-t table_file] table of compositions to be substituted into melts_file \n".
	      " [-o output_file] filename for program output \n".
	      " [-p path] path of directory in which to put output files\n".
	      "              (if not the same as working directory) \n".
	      " [-c column_file] run 'column_pick.command < column_file' before moving files\n".
	      " [-l log_file] file to record environment and initial conditions \n".
	      "              (default logfile.txt) \n".
	      " P, T and fO2 may be initialised in the melts_file, settings_file or program:\n".
	      "  - only settings_file P, T, fO2 will be recorded in the log_file.\n".
	      "  - settings_file P, T, fO2 will overwrite melts_file values.\n\n".
	      " This is run_alphamelts.command 1.9; use it with alphamelts 1.9 or updates 1.9.X.\n\n");
	next;
    }

    $path = $batch = '';
    $run_path = $0; # $0 is the name (and path) of the run_alphamelts script
#    $run_path =~ s/run_alphamelts.command//;
    $run_path =~ s/run_alphameltsExoC.command//;
    $run = (-f "$run_path/with-readline") ? 'with-readline ' : '';
    $program = 'alphamelts';
    
    $log_file = 'logfile.txt';
    $table_file = $melts_file = $batch_file = $column_file = '';
    $out_file = 'alphaMELTS_tbl.txt';

    $rs = $/; # expected record separator i.e. \r\n or \n
    $delim = '/';
    $windows = '';
    if (exists $ENV{'OS'}) {
	if ($ENV{'OS'} =~ /windows/i) {
	    $windows = 1;
	}
    }
    $delim = '\\' if ($windows);

    if (grep /^-/, @argv2) {

	my $temp;

	# check command line switches
	map {

	    if (/^-/) {

		$temp = $_;

	    }
	    else {

		s/^\"//; # Trim any leftover quotes from start and end of list
		s/\"$//;
		# Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
		s/\"\\\"\"/\'/g || s/\\\"/\'/g || s/\"/\'/g;
		
		s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')
		$com{$temp} = $_;

	    }

	} @argv2;

	$path = $com{'-p'} if (exists $com{'-p'});
	$path .= $delim if ($path); # path and / for good measure

	$run = 'gdb ' if (grep /^-g/, @argv2);
	$run = 'lldb ' if (grep /^-d/, @argv2);
	$run = 'valgrind --leak-check=full --track-origins=yes ' if (grep /^-v/, @argv2);

	$batch = grep /^-(a|b)/, @argv2;
	$batch_file = ($com{'-b'})  if (exists $com{'-b'});

	$in_file = ($com{'-f'})  if (exists $com{'-f'});
	$melts_file = ($com{'-m'})  if (exists $com{'-m'});
	$table_file = ($com{'-t'})  if (exists $com{'-t'});
	$out_file = ($com{'-o'})  if (exists $com{'-o'});
	$log_file = ($com{'-l'})  if (exists $com{'-l'});
	$column_file = ($com{'-c'})  if (exists $com{'-c'});
	
	# check command line switches
	if ((grep /^-[^hpfbamtolcgdv]/, @argv2) || exists $com{'-g'} || exists $com{'-d'} || exists $com{'-v'}
	    || exists $com{'-a'}) {
	    grep {warn "RUN_ALPHAMELTS.COMMAND ERROR: Unknown option \'$_\'\n" if $_ !~ /^-[hpfbamtolcgv]/} @argv2;
	    grep {warn "RUN_ALPHAMELTS.COMMAND ERROR: Unknown option \'$_ $com{$_}\'\n" if exists $com{$_}} ('-g', '-d', '-v', '-a');
	    next;
	}

    }
    elsif (@argv2) {
	grep {warn "RUN_ALPHAMELTS.COMMAND ERROR: Unknown option \'$_\'\n"} @argv2;
	next;
    }
    else {
	unshift @argv2, "\(none\)";
    }

    (-d $path || mkdir "$path") if ($path);

    $log_file = "$path$log_file";
    unless (open (LOGFILE, ">$log_file")) {
	warn "RUN_ALPHAMELTS.COMMAND WARNING: Cannot open log_file \"$log_file\" ($!)\n";
	next;
    }
    unless (open (INPUT, "<$in_file")) {
	warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open settings_file \"$in_file\" ($!)\n";
	next;
    }

    @oldlines = <INPUT>; # read whole file into array of lines - O.K. as it's not very big
    @oldlines = split /[\r\n]+/, (join '', @oldlines); # line endings may not be correct

    close(INPUT);

    print LOGFILE map {$_ .= "\n"} @oldlines;
    map {s/^\s*//; s/\=//; s/!.*//; chomp} @oldlines; # remove leading white space and trailing comments
    @oldlines =  grep /.+/, @oldlines;

    $subsol =  grep {
	if  (/Subsolidus Phases/i) {
	    @phases = split /\s+/, (split /\:\s+/)[1];  # add carriage return to each phase
	    map {$_ .= "\n"} @phases;  # strip out white space first
	    push @phases, "x\n"; # x to exit phase list in alphamelts
	}
    } @oldlines;
    
    if($melts_file) {

	unless (open (OLDMELTS, "<$melts_file")) {
	    warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open melts_file \"$melts_file\" ($!)\n";
	    next;
	}

	my @newlines = <OLDMELTS>; # read whole melts file into array of lines - O.K. as it's not very big	
	@newlines = split /[\r\n]+/, (join '', @newlines); # line endings may not be correct
	map {$_ .= "\n"} @newlines;

	close OLDMELTS;
	unless ((rename $melts_file, "$melts_file\_bak") && (open (NEWMELTS, ">$melts_file"))) {
	    warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open melts_file \"$melts_file\" ($!)\n";
	    next;
	}

	# ignore environment variables and certain other commands for now
	@lines = grep !(/(ALPHAMELTS|Subsolidus)/i), @oldlines;

	map {
	    $line = "$_\n"; # remaining lines should be initial P, T, FO2 or mass
	    print "$line"; # write line to screen
	    if (/Suppress:/i) {
		if(/none/i) {
		    @newlines = grep !(/Suppress/i), @newlines; # no phase names contain 'none'
		} 
		else {
		    # need to think about numbers for feldspar / pyx !?
		    push @newlines, $line unless grep /$line/, @newlines;
		}
	    }
	    else {
		# substitute old line in melts file for new line from input file
		# or add new line if not previously set in melts file
		push @newlines, $line unless grep {s/$_/$line/i if $_ =~ (split /\:/, $line)[0]} @newlines;
	    }
	} @lines;
	    
	print NEWMELTS @newlines; # write new melts file
	close NEWMELTS;
	
	if($table_file) {

	    unless (open (TABLE, "<$table_file")) {
		warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open table_file \"$table_file\" ($!)\n";
		next;
	    }

	    @lines = <TABLE>; # read whole file into array of lines - O.K. if it's not very big
	    @lines = grep /.+/, (split /[\r\n]+/, (join '', @lines)); # line endings may not be correct

	    close TABLE;

	    @used_oxides = (); $table_row = 1; # anything non-zero would work
	    map {
		
		s/^\s*//;
		if(@used_oxides) {

		    my (%values, $label, $value);

		    my $new_file = $melts_file;
		    $table_row++;
		    $new_file =~ s/\./$table_row\./ || ($new_file .= $table_row);

		    open (NEWMELTS, ">$new_file");

		    grep {# add row number to title
			s/\($table_file row .*/\($table_file row $table_row\)/i || 
			    s/^Title.*/$& \($table_file row $table_row\)/i
		    } @newlines;

		    $label = 0;
		    map {$values{$used_oxides[$label++]} = $_} split /[ \t,]+/;

		    foreach $label (@used_oxides) {
			$value = $values{$label};
			grep {# initial composition or trace (check for equality so P does not match Pb etc.)...
			    /^(Initial|Log)/i && ((lc $label eq lc ((split)[2])) ? s/\s$label.*/ $label\ $value/i
						  : s/\s$label:.*/ $label:\ $value/i) # ... or any other initial value
			} @newlines;
		    }
		    
		    print NEWMELTS @newlines; # write new melts file
		    close NEWMELTS;

		}
		else {

		    my ($label, @test);
		    @used_oxides = split /[ \t,]+/;

		    foreach $label (@used_oxides) { #at least one match
			@test = grep {/^(Initial|Log)/i && ((lc $label eq lc ((split)[2])) || /$label:\ .+/i)} @newlines;
			if(@test) {
			    $table_row = 0; #use as a counter of lines after header
			    next;
			}
		    }
		    @used_oxides = () if ($table_row); # not yet got to header

		}
	    } @lines;
	    
	}

	# emulate alphaMELTS 2 behaviour
	my (@temp, $init, $final, $inc);
	$init = ''; $final = ''; $inc = '';
	map {
	    if (/temperature/i) {
		chomp;
		my $val = (split /\:\s+/)[1];
		$init = $val if (/initial/i);
		$final = $val if (/final/i);
		$inc = $val if (/increment/i);
	    }
	} @newlines;
	if ($init && ($final > $init)) {
	    $ENV{"ALPHAMELTS_MAXT"} = $final;
	    $ENV{"ALPHAMELTS_DELTAT"} = abs $inc if ($inc);
	}
	elsif ($final && ($init > $final)) {
	    $ENV{"ALPHAMELTS_MINT"} = $final;
	    $ENV{"ALPHAMELTS_DELTAT"} = -1 * abs $inc if ($inc);
	}
	

	$init = ''; $final = ''; $inc = '';
	map {
	    if (/pressure/i) {
		chomp;
		my $val = (split /\:\s+/)[1];
		$init = $val if (/initial/i);
		$final = $val if (/final/i);
		$inc = $val if (/increment/i);
	    }
	} @newlines;
	if ($init && ($final > $init)) {
	    $ENV{"ALPHAMELTS_MAXP"} = $final;
	    $ENV{"ALPHAMELTS_DELTAP"} = abs $inc if ($inc);
	}
	elsif ($final && ($init > $final)) {
	    $ENV{"ALPHAMELTS_MINP"} = $final;
	    $ENV{"ALPHAMELTS_DELTAP"} = -1 * abs $inc if ($inc);
	}
		    
	map {	    
	    if (/Mode\:/i) {
		$ENV{"ALPHAMELTS_FRACTIONATE_SOLIDS"} = 1 if (/fractionate solids/i);
		$ENV{"ALPHAMELTS_FRACTIONATE_WATER"} = 1 if (/fractionate fluids/i);
		if (/fractionate liquids/) {
		    $ENV{"ALPHAMELTS_CONTINUOUS_MELTING"} = 1;
		    $ENV{"ALPHAMELTS_FRACTIONATE_SECOND_LIQUID"} = 1;
		}
		$ENV{"ALPHAMELTS_MULTIPLE_LIQUIDS"} = 1 if (/multiple liquids/i);
		
		$ENV{"ALPHAMELTS_MODE"} = "isenthalpic" if (/isenthalpic/i);
		$ENV{"ALPHAMELTS_MODE"} = "isentropic" if (/isentropic/i);
		$ENV{"ALPHAMELTS_MODE"} = "isochoric" if (/isochoric/i);
	    }		
	} @newlines;

	@newlines = ();
    }
    elsif ($table_file) {
	warn "RUN_ALPHAMELTS.COMMAND ERROR: Please provide a melts_file to use as template for the table_file\n";
	next;
    }
    elsif ($batch && !$batch_file) {
	warn "RUN_ALPHAMELTS.COMMAND ERROR: Please provide a melts_file to use in automatic mode\n";
	next;
    }
    @lines = ();

    # Make sure settings file overrides .melts
    map {

	if (/^ALPHAMELTS/) {
	    %ENV = (
		%ENV, # existing environment variables
		split # split by white space into key and value
		) 
	}
	elsif (/^ADIABAT/i) {
	  warn "'ADIABAT' environment variables are no longer recognized.\n";
	  warn "Please replace 'ADIABAT' with 'ALPHAMELTS'.\n";
	  next;
	}

    } @oldlines;
    @oldlines = ();
    
    if (exists $ENV{'ALPHAMELTS_INTEGRATE_FILE'}) {
	my $temp_file = $ENV{'ALPHAMELTS_INTEGRATE_FILE'};
	$ENV{'ALPHAMELTS_INTEGRATE_FILE'} = "$path$temp_file";
    }
    
    map {print "$_ $ENV{$_}\n" if /^ALPHAMELTS/;} keys %ENV;

# ABOUT TO RUN ALPHAMELTS_1PH ******************

# these commands could be used to modify this file to run other batch procedures

    %com = (
	quit =>  0,
	read =>  1,
	twiddle => 2,
	single => 3,
	exe =>  4,
	super => 1,
	subs =>  0,
	);

    if($batch) {

	if(!$batch_file) {

	    my ($i);
	    $batch_file = "auto_batch.txt";

	    open (BATCH, ">$batch_file");

	    print BATCH "run\n" if ($run && ($run !~ /with-readline/));

	    if ($table_file && $table_row) {
		for ($i = 1; $i <= $table_row; $i++) {
		    my $new_file = $melts_file;
		    $new_file =~ s/\./$i\./ || ($new_file .= $i);
		    
		    print BATCH "$com{'read'}\n\"$new_file\"\n";
		    print BATCH "$com{'single'}\n";

		    print BATCH ($subsol) ? "$com{'subs'}\n" : "$com{'super'}\n";
		    print BATCH @phases if($subsol);

		    print BATCH "$com{'exe'}\n";
		}
	    }
	    else {
		print BATCH "$com{'read'}\n\"$melts_file\"\n";
		print BATCH "$com{'single'}\n";

		print BATCH ($subsol) ? "$com{'subs'}\n" : "$com{'super'}\n";
		print BATCH @phases if($subsol);

		print BATCH "$com{'exe'}\n";
	    }

	    print BATCH "$com{'quit'}\n";
	    close BATCH;

	}

	# $program does not have spaces etc. but $run_path might
	#((-f "$run_path$program") && !(system "$run\"$run_path$program\" < $batch_file")) || 
	#    (!(-f "$run_path$program") && !(system "$run$program < $batch_file")) ||
	#    warn "RUN_ALPHAMELTS.COMMAND WARNING: alphamelts may have crashed!\n";

	# $program does not have spaces etc. but $run_path might
	((-f "/Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/alphamelts_macosx64") && !(system "/Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/alphamelts_macosx64 < /Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/ExoC/ExoCbatchLith.txt")) || 
	    (!(-f "$run_path$program") && !(system "$run$program < $batch_file")) ||
	    warn "RUN_ALPHAMELTS.COMMAND WARNING: alphamelts may have crashed!\n";

    }
    else {

	# print relevant file names before execution if using alphamelts interactively
	if ($table_file) {
	    my $new_file = $melts_file;
	    $new_file =~ s/\./1\./ || ($new_file .= '1');
	    print "\nPlease use melts files: $new_file, ";
	    $new_file =~ s/1\./2\./ || s/1$/2/;
	    print "$new_file etc.\n\n";
	}
	else {
	    ($melts_file) ? print "\nPlease use melts file: $melts_file\n\n" : print "\n\n";
	}

	# Try $run_path first as can test for the file; then try path (e.g. for Mac double-click)
	#((-f "$run_path$program") && !(system "$run\"$run_path$program\"")) ||
	#    (!(-f "$run_path$program") && !(system "$run$program")) ||
	#    warn "RUN_ALPHAMELTS.COMMAND WARNING: alphamelts may have crashed!\n";

	# Try $run_path first as can test for the file; then try path (e.g. for Mac double-click)
	((-f "/Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/alphamelts_macosx64") && !(system "/Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/alphamelts_macosx64")) ||
	    (!(-f "$run_path$program") && !(system "$run$program")) ||
	    warn "RUN_ALPHAMELTS.COMMAND WARNING: alphamelts may have crashed!\n";
	
    }

# FINISHED RUNNING ALPHAMELTS_1PH **********************

    print "\n\n";

    {

	my ($me, @months, $year, $month, $date, $hour, $min);

	($me = $ENV{'USERNAME'}) || ($me = $ENV{'USER'}); # system environment variable
	chomp $me;

	@months = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec');

	($min, $hour, $date, $month, $year) = (localtime)[1..5];

	$year += 1900;
	$month = $months[$month];
	$hour = sprintf "%2.2d", $hour;
	$min = sprintf "%2.2d", $min;

	# record run conditions in log file
	print LOGFILE "\n\! alphamelts run by $me \n".
	    "! Time: $month $date $year at $hour\:$min \n".
	    "! Command line: @argv2 \n".
	    "! Settings file: $in_file\n";
	close LOGFILE;

    }

######################################################################

    # Test solids first so don't waste time processing old file
    unless (-f "Phase_main_tbl.txt") {
	warn "RUN_ALPHAMELTS.COMMAND WARNING: Cannot find output file \"Phase_main_tbl.txt\" ($!)\n";
	warn "Please check that alphamelts ran properly.\n";
	next;
    }

    if (open (INPUT, "<Phase_main_tbl.txt")) {
	@lines = <INPUT>;
	close (INPUT);
	unless (grep /(Pressure\s+)(.+)(Temperature.*)/, @lines) {
	    # probably means the Solids file has already been processed
	    warn "RUN_ALPHAMELTS.COMMAND WARNING: Incorrect format in output file \"Phase_main_tbl.txt\".\n";
	    warn "File has already been processed? Please check that alphamelts ran properly.\n";
	    next;
	}
    }
    else {
	warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open out_file \"$out_file\" ($!)\n";
	next;
    }

    my @out_files = ("System_main_tbl.txt",
		     "Liquid_comp_tbl.txt",
		     "Phase_main_tbl.txt",
		     "Phase_mass_tbl.txt",
		     "Phase_vol_tbl.txt",
		     "Solid_comp_tbl.txt",
		     "Bulk_comp_tbl.txt");

    unless (open (OUTPUT, ">$out_file")) {
	warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open out_file \"$out_file\" ($!)\n";
	next;
    }

    foreach $out_file (@out_files) {

	# now start processing output; first check file exists
	unless ((-f $out_file) && (rename $out_file, "$out_file\_bak") && 
		(open (INPUT, "<$out_file\_bak")) && (open (TABLE, "+>$out_file"))) { 
	    warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open out_file \"$out_file\" ($!)\n";
	    next;
	}

	if (lc $out_file eq lc "Phase_main_tbl.txt") {

	    my (@columns, $label, @string, $pt, $key, $val, %therm_comp, $el_list);
    
	    $title = '';
	    while (<INPUT>) {

		chomp;
		if ($title) {
		    print TABLE "$title\n" if ($title); # print the title
		    $title = ''; 
		}
		elsif (/^Title\:/) {
		    $title = $_;  # print the title next time (above)
		}
		elsif (/^Pressure/) { # get the new P and T value

		    $pt = $_;
		    $pt =~ s/(Pressure\s+)(.+)(Temperature.*)/$2/;
		    s/(Pressure\s+)(.+)(Temperature\s+)//;
		    @columns = split /\s+/; # open once to get element list
		    $pt .= shift @columns;
		    $el_list = "@columns";

		}
		else { 
		
		    s/\ oxide/_oxide/g;
		    s/\ ss/_ss/g;
		    @columns = split /\s+/;
		    $label = shift @columns;
		    
		    if($pt) {
			# therm_comp is a hash, with the phase name as the key and an array of lines as the corresponding entry
			push @{$therm_comp{"$label"}}, "$pt @columns\n"; # add a line to the appropriate array
		    }
		    else {
			# probably means the Solids file has already been processed
			warn "RUN_ALPHAMELTS.COMMAND WARNING: Incorrect format in output file \"$out_file\".\n";
			warn "File has already been processed? Please check that alphamelts ran properly.\n";
			next;
		    }

		}
	    }
    
	    $pt = "Pressure Temperature mass S H V Cp";

	    while (($key, $val) = each %therm_comp) {
		if ($key =~ /^liquid/) { # liquid(s) first
		    @string = @$val;
		    print TABLE "\n$key thermodynamic data and composition:\n$pt viscosity $el_list Mg#\n";
		    print TABLE @string;
		}
	    }
	    while (($key, $val) = each %therm_comp) {
		unless ($key =~ /^liquid/) { # solids after
		    @string = @$val;
		    if ($key =~ /(pyroxene|amphibole)/) {
			print TABLE "\n$key thermodynamic data and composition:\n$pt structure formula $el_list\n";
		    }
		    elsif ($key =~ /\_[0-9]/) { # solid solution
			print TABLE "\n$key thermodynamic data and composition:\n$pt formula $el_list\n";
		    }
		    else { # pure phase
			print TABLE "\n$key thermodynamic data and composition:\n$pt formula\n";
		    }
		    print TABLE @string;
		}
	    }

	}
	elsif (lc $out_file eq lc "Phase_mass_tbl.txt") {
	    while (<INPUT>) {
		if (/^Pressure/) {
		    s/\ oxide/_oxide/g;
		    s/\ ss/_ss/g;
		}
		print TABLE;
	    }
	}
	else {
	    while (<INPUT>) {
		print TABLE;
	    }
	}

	# slurp
	seek (TABLE, 0, 0);
	@lines = <TABLE>;
	print OUTPUT @lines;
	print OUTPUT "\n";

	close (INPUT);
	close (TABLE);
	
	unless (unlink "$out_file\_bak")  {
	    warn "Could not delete backup file \"$out_file\_bak\" ($!)\n";
	    next;
	}

    }

    # now start processing the trace file if appropriate; check it exists and back up
    if (exists $ENV{'ALPHAMELTS_DO_TRACE'} || exists $ENV{'ALPHAMELTS_DO_TRACE_H2O'}) {

	my (@columns, $label, @string, $pt, $key, $val, %therm_comp, $el_list);
    
	$out_file = "Trace_main_tbl.txt";
	unless ((-f $out_file) && (rename $out_file, "$out_file\_bak") && 
		(open (INPUT, "<$out_file\_bak")) && (open (TABLE, "+>$out_file"))) {
	    warn "RUN_ALPHAMELTS.COMMAND ERROR: Cannot open out_file \"$out_file\" ($!)\n";
	    next;
	}

	$title = '';
	while (<INPUT>) {

	    chomp;
	    if ($title) {
		print TABLE "$title\n" if ($title); # print the title
		$title = ''; 
	    }
	    elsif (/^Title\:/) {
		$title = $_;  # print the title next time (above)
	    }
	    elsif (/^Pressure/) { # get the new P and T value

		$pt = $_;
		$pt =~ s/(Pressure\s+)(.+)(Temperature.*)/$2/;
		s/(Pressure\s+)(.+)(Temperature\s+)//;
		@columns = split /\s+/; # open once to get element list
		$pt .= shift @columns;
		$el_list = "@columns";

	    }
	    else { 
	    
		s/\ oxide/_oxide/g;
		s/\ ss/_ss/g;
		@columns = split /\s+/;
		$label = shift @columns;
		
		if($pt) {
		  # therm_comp is a hash, with the Partition, Bulk or phase name as the key and an array of lines as the corresponding entry
		  push @{$therm_comp{"$label"}}, "$pt @columns\n"; # add a line to the appropriate array
		}
		else {
		  # probably means the Traces file has already been processed
		  warn "RUN_ALPHAMELTS.COMMAND WARNING: Incorrect format in output file \"$out_file\".\n";
		  warn "File has already been processed? Please check that alphamelts ran properly.\n";
		  next;
		}

	    }
	}
    
	$pt = "Pressure Temperature mass ";

	$val = $therm_comp{'Partition'}; # print out in same order as before
	@string = @$val;
	print TABLE "\nPartition Coefficients:\n$pt$el_list\n";
	print TABLE @string;
	$val = $therm_comp{'Bulk'};
	@string = @$val;
	print TABLE "\nBulk Trace:\n$pt$el_list\n";
	print TABLE @string;
	$val = $therm_comp{'Solid'};
	@string = @$val;
	print TABLE "\nSolid Trace:\n$pt$el_list\n";
	print TABLE @string;
	if (exists $therm_comp{'Liquid'}) {
	    $val = $therm_comp{'Liquid'};
	    @string = @$val;
	    print TABLE "\nLiquid Trace:\n$pt$el_list\n";
	    print TABLE @string;
	}
	while (($key, $val) = each %therm_comp) {
	    unless ($key =~ /Partition|Bulk|Solid|Liquid/) { # the rest
		@string = @$val;
		print TABLE "\n$key trace:\n$pt$el_list\n";
		print TABLE @string;
	    }
	}
	
	# slurp
	seek (TABLE, 0, 0);
	@lines = <TABLE>;
	print OUTPUT @lines;
	print OUTPUT "\n";

	close(INPUT);
	close(TABLE);

	unless (unlink "$out_file\_bak") {
	    warn "Could not delete backup file \"$out_file\_bak\" ($!)\n";
	    next;
	}

    }
    close(OUTPUT);

}
continue {

    # If forgot -c or -p switches don't have to rerun whole program - just choose menu option '0'
    if ($column_file) {
	$program = 'column_pick.command';
	((-f "$run_path$program") && !(system "$run\"$run_path$program\" < $column_file")) || 
	    (!(-f "$run_path$program") && !(system "$run$program < $column_file")) ||
	    warn "RUN_ALPHAMELTS.COMMAND WARNING: column_pick may not have run properly!\n";
    }
    if ($path) {
	$program = ($windows) ? 'move /Y' : 'mv -f';
	my @files = glob '*_tbl.txt';
	for my $file (@files) {
	    system "$program $file \"$path\""; 
	}
    }

    unless (@ARGV) {
	%com = ();
	print ("Run again (y or n)?\n");
	$_ = <STDIN>;
	chomp;
	# Mostly stays open if user accidently drags and drops before answering...
	@argv2 = () unless (/^$/ || /^n/i); # shuts if no answer or begins with n/N.
    }

}
