#!/usr/bin/perl -w

use strict;

my ($savefile, @argv2, $path, $delim, $header, $file, $table, $tablename, @pt1, @pt2, @pt3);
my (@labels, @columns, %values, @output, $lineno, @lines, $windows, $infile, $outfile);

$| = 1; # print buffer immediately
$delim = ' '; # default delimiter is space
$header = 'excel'; # default format is two header lines

$file = '';
$table = '';
$tablename = '';
@pt1 = @pt2 = @pt3 = ();

$lineno = 2;
@lines = ();

$windows = 0;
if(exists($ENV{'OS'})) {
  if($ENV{'OS'} =~ /Windows/i) {
    $windows = 1;
  }
}

$infile = '';
$outfile = '';
until ($outfile && !$infile) {

    if (@ARGV) {

 	@argv2 = @ARGV;
	$_ = $argv2[$#argv2];
	if (/\>./) {
	    $outfile = '>'; # Redirect to Matlab
	    pop @argv2;
	}
	else {
	    $outfile = @ARGV;
	}

    }
    else {

	if ($infile) {
	    $_ = <INFILE>;
	}
	else {
	    print ("Usage: column_pick.command column_list_files > table_file\n");
	    print ("Please enter column_list_files, separated by spaces, and / or\n");
	    print ("type '>' and the table_file name to finish and write output.\n");
	    $_ = <STDIN>;
	}
	
	chomp;

	unless (/.+/) { # Reached empty line or end of file
	    $infile = '';
	    next;
	}

	if (/\>/) {
	    $outfile = $_;
	    $outfile =~ s/^(.*)(\>)(\s*)(.+)$/$4/;
	    $outfile =~ s/\s+$//;
	}
	s/^(.*)(\>)(.*)$/$1/; # Trim output file name and any trailing white space

	if (/\</) {

	    $infile = $_;
	    $infile =~ s/^(.*)(\<)(\s*)(.+)$/$4/;
	    $infile =~ s/\s+$//;
	    $infile =~ s/^(\"|\')//; # Trim any leftover quotes from start and end of file name
	    $infile =~ s/(\"|\')$//;
	    # Fix any apostrophes in the name (if had single quotes or no quotes originally)
	    $infile =~ s/\'\\\'\'/\'/g || $outfile =~ s/\\\'/\'/g;
	    $infile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

	    unless (open(INFILE, "<$infile")) {
		warn "Error: could not open file \"$infile\" : $!\n";
		$infile = '';
		next;
	    }

	}

	s/\s+$//;
	s/\'/\"/g; # Windows can only deal with name enclosed in double quotes; *nix can take single or double
	# For a list of files inside quotes it follows that the white space is inside quotes (e.g. '" "' in '"file 1" "file 2"')
	# Otherwise split by spaces and adjust later if any are escaped (e.g. ' ' versus '\ ' in 'file\ 1 file\ 2')
	@argv2 = (/[^\\]\"/) ? split /\"\s+\"/ : split /\s+/; # ... but don't match escaped apostrophes (currently '\"')

    }

    FOREACH: foreach my $colfile (@argv2) {

   	$colfile =~ s/^\"//; # Trim any leftover quotes from start and end of list
	$colfile =~ s/\"$//;
	# Fix any apostrophes in the name (if had single quotes or no quotes or double quotes originally)
 	$colfile =~ s/\"\\\"\"/\'/g || $colfile =~ s/\\\"/\'/g || $colfile =~ s/\"/\'/g;

        # A trailing backslash indicates an escaped space (except on Windows where it means a folder, not a file)
	$colfile = $savefile.' '.$colfile if ($savefile);
        $savefile = ($colfile =~ s/\\$//) ? $colfile : '';
        next if ($savefile);

	$colfile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

	unless (open(COLUMNLIST, "<$colfile") || ($outfile eq '>')) {
	    warn "Error: could not open file \"$colfile\" : $!\n";
	    last FOREACH;
	}
	($outfile eq '>') || warn "*** \"$colfile\" ***\n\n";

	$path = '';

      MAIN: while(<COLUMNLIST>) {

	  # Process the column_list_file triplets
	  s/!.*//;
	  chomp;
	  if(/^Name/) {

	      $path = $colfile;
	      $colfile = $_;
	      $colfile =~ s/Name(:*)(\s+)//;
	      $colfile =~ s/^[\'\"]//; # Trim any quotes from start and end
	      $colfile =~ s/[\'\"]$//;
	      ($outfile eq '>') || warn ("Filename set to: \"$colfile\"\n");

	  }
	  elsif(/^Delimiter/i) {

	      $delim = (split)[1];
	      $delim =~ s/[\'\"]//g; # strip quotes
	      if(($delim =~ /\\t/) || (lc $delim eq 'tab')) {
		  $delim = "\t";
		  ($outfile eq '>') || warn ("Delimiter set to: '\\t'\n\n");
	      }
	      else {
		  if(lc $delim eq 'space') {
		      $delim = ' ';
		  }
		  elsif(lc $delim eq 'comma') {
		      $delim = ',' ;
		  }
		  ($outfile eq '>') || warn ("Delimiter set to: '$delim'\n\n");
	      }
	      
	  }
	  elsif(/^Header/i) {

	      $header = lc ((split)[1]);
	      $header =~ s/[\'\"]//g; # strip quotes
	      ($outfile eq '>') || warn ("Header type set to: '$header'\n\n");
	      
	  }
	  elsif(/^File/i) {
	      
	      s/^File(s*)(\:*)(\s+)//; # The alphaMELTS (or other) output file
	      s/^[\'\"]//; # Trim any quotes from start and end
	      s/[\'\"]$//;
	      $file = $_;

	      if ($path) {
		  $_ = $path;
		  s/$colfile/$file/;
		  while (/\.\./) {
		      ($windows) ? (s/\\[^\.]+\\\.\.\\/\\/ || s/[^\.]+\\\.\.\\//) : 
			  (s/\/[^\.]+\/\.\.\//\// || s/[^\.]+\/\.\.\///);
		  }
		  $file = $_;
	      }
	      ($outfile eq '>') || warn ("File set to: \"$file\"\n");

	  }
	  elsif(/^Table/i) {

	      s/[a-z\:\s]+[\'\"]//i; # Strip up to first quote, if present
	      $table = (/[\'\"]/) ? $_ : (split)[1]; # First word(s) of the table name in that file
	      $table =~ s/[\'\"]//g; # strip the other quote
	      ($outfile eq '>') || warn ("Table name set to: '$table'\n");

	  }
	  elsif(/^Variable/i) {

	      s/[\'\"]//g; # strip quotes
	      @pt1 = @pt2 = @pt3 = (); # empty arrays to unset the variable values
	      @columns = split; # list of variable names
	      shift @columns; # remove the 'Variables:' part

	      if(@columns) {

		  # Process text file to pick up variable values for comparison later ...
		  unless (open(FILE, "<$file") || ($outfile eq '>')) {
		      warn ("Error: could not open file \"$file\" at line $. : $!\n");
		      last FOREACH;
		  }
		  ($outfile eq '>') || warn "... now reading \"$file\" to look for variables\n";
		  
		VARIABLES: while(<FILE>) {

		    chomp;
		    if (!$tablename) {

			# Skip forward to the table whose title begins with the given word(s),
			# with or without trailing underscore zero for first appearance of phases.
			next unless (/^$table\_0/i || ("$_\_0" =~ /^$table/i) || /^$table/i);
			$tablename = $_;

		    }
		    elsif (/\s[0-9-]/) {

			# For alphaMELTS output, no column names begin with a number or minus sign,
			# whereas a line of data will always have at least one number or '---'
			my $i = 0;
			map {$values{$labels[$i++]} = $_} split /\s+/; # one value for each matched column label

			# Append the data to array(s) of variable value(s)
			# Arrays are set up below and zeroth position is column label
			push @pt1, $values{$pt1[0]} if (@pt1); 
			push @pt2, $values{$pt2[0]} if (@pt2);
			push @pt3, $values{$pt3[0]} if (@pt3);

		    }
		    elsif(/.+/) {

			# The line is not empty so, by elimination, it must be the column names
			@labels = split /\s+/;

			# Check that the named variables are actually present (doesn't matter where)
			# Case and '_0' do not matter but otherwise it must be an exact match (i.e. 'pressure' not 'P')
			foreach my $column (@columns) {
			    my $label = (grep {(lc $column eq lc $_) || (lc "$column\_0" eq lc $_) ||
						   (lc $column eq lc "$_\_0")} @labels)[0];
			    
			    if ($label) {
				$values{$column} = $label;
			    }
			    else {
				($outfile eq '>') ||  
				    warn ("Warning: could not find variable '$column' in table '$table', file \"$file\" at line $.!\n");
				$values{$column} = '';
			    }

			}

			# Set up array(s) for the variable(s) and put column label in first place
			push @pt1, $values{$columns[0]};
			push @pt2, $values{$columns[1]} if (@columns > 1);
			push @pt3, $values{$columns[2]} if (@columns > 2);

			# This space will be used for the column name in the matched table (i.e. 'mass' if not the same).
			push @pt1, '' if (@pt1);
			push @pt2, '' if (@pt2);
			push @pt3, '' if (@pt3);
	    
			foreach my $column (@columns) {
			    $column = $values{$column}; # list of successfully matched column labels
			}
			($outfile eq '>') ||  warn ("Variables set to: @columns\n\n");

		    }
		    else{
			# blank line signifies end of relevant table so jump out of loop
			last VARIABLES;
		    }
		    %values = ();
		}

		  %values = (); # hit end of file, either successfully or not
		  warn ("Warning: could not find table '$table' in file \"$file\" at line $.!\n") unless (($outfile eq '>') || $tablename);
		  $lineno = 2;
		  close FILE;
		  $tablename = ''; # reset for next go

	      }
	      else {

		  # ... or if the 'Variables:' line is empty just leave the variables unset and continue
		  ($outfile eq '>') || warn ("Variables are now unset ...\n\n");

	      }


	  }
	  elsif(/^Column/i) {

	      s/[\'\"]//g; # strip quotes
	      @columns = split; # list of column names
	      shift @columns; # remove the 'Variables:' part

	      # Process text file to pick up matching data
	      unless (open(FILE, "<$file") || ($outfile eq '>')){ 
		  warn ("Error: could not open file \"$file\" at line $. : $!\n");
		  last FOREACH;
	      }
	      ($outfile eq '>') || warn "... now reading \"$file\" to look for columns\n"; 

	    COLUMNS: while(<FILE>) {

		chomp;
		if (!$tablename) {

		    # Skip forward to the table whose title begins with the given word(s),
		    # with or without trailing underscore zero for first appearance of phases.
		    next unless (/^$table\_0/i || ("$_\_0" =~ /^$table/i) || /^$table/i);

		    my $i = @columns;
		    if ($header eq 'matlab') {
			$_ = $table;
			$tablename = (split)[0]; # use just first word of table name in printed columnn names
		    }
		    else {
			$tablename = (/^$table\_0/i) ? $table : $_; # for output print '_0' only if given
			# encase full table name in quotes so it will be imported into a single spreadsheet cell
			$lines[0] .= "\"$tablename\"";
		    }
		    $lines[0] .= $delim x $i; # pad so that next table name will be correctly aligned

		}
		elsif (/\s[0-9-]/) {

		    # For alphaMELTS output, no column names begin with a number or minus sign,
		    # whereas a line of data will always have at least one number or '---'
		    my $line = $_;
		    my $i = 0;
		    map {$values{$labels[$i++]} = $_} split /\s+/; # one value for each matched column label
		    $values{$labels[$i-1]} = '' if ($line =~ /---/); # put in empty value instead of '---'

		    # For debugging purposes:
		    #my $warning = "\n";
		    #map {$warning .= "$_ $values{$_} "} keys %values;
		    #$warning .= "\n".((@pt1) ? "$pt1[1] $pt1[$lineno]" : '').' '.((@pt2) ? "$pt2[1] $pt2[$lineno]" : '').' '
		    #  .((@pt3) ? "$pt3[1] $pt3[$lineno]" : '');
		    #warn $warning;

		    # if the @pt1 array is empty then '!@pt1' is 'true'; otherwise make sure the variable value matches etc.
		    my $match = '';
		    until ($match) { # evaluate loop at least once

			$match = ((!@pt1 || ($values{$pt1[1]} eq $pt1[$lineno])) &&
				  (!@pt2 || ($values{$pt2[1]} eq $pt2[$lineno])) &&
				  (!@pt3 || ($values{$pt3[1]} eq $pt3[$lineno])));
			
			if($match) {
			    foreach my $column (@columns) { # all columns, including the missing ones
				push @output, $values{$column};
				$values{$column} = ''; # reinitialise (will stay empty for missing columns and if have '---')
			    }
			    $lines[$lineno++] .= (join $delim, @output).$delim;
			}
			else {
			    # pad those lines for which there are variables set but no match
			    my $i = @columns; # number of columns
			    $lines[$lineno++] .= $delim x $i;
			}

		    }

		}
		elsif(/.+/) {

		    # The line is not empty so, by elimination, it must be the column names
		    @labels = split /\s+/;

		    # If the @pt1 array is empty then '!@pt1' is 'true'; otherwise make sure the variable name matches etc.
		    # Put 'mass' in column name if necessary
		    !@pt1 || ($pt1[1] = ((lc $pt1[0] eq lc $table) || (lc "$pt1[0]\_0" eq lc $table) || (lc $pt1[0] eq lc "$table\_0")) ?
			      'mass' : $pt1[0]);
		    !@pt2 || ($pt2[1] = ((lc $pt2[0] eq lc $table) || (lc "$pt2[0]\_0" eq lc $table) || (lc $pt2[0] eq lc "$table\_0")) ?
			      'mass' : $pt2[0]);
		    !@pt3 || ($pt3[1] = ((lc $pt3[0] eq lc $table) || (lc "$pt3[0]\_0" eq lc $table) || (lc $pt3[0] eq lc "$table\_0")) ?
			      'mass' : $pt3[0]);
  
		    # Print header line before substituting true column name for given label / format
		    if ($header eq 'matlab') { # put (short) table name
			$lines[1] .= (join "_$tablename$delim", @columns)."_$tablename$delim";
		    }
		    else { # table name already printed
			$lines[1] .= (join $delim, @columns).$delim;
		    }

		    foreach my $column (@columns) {
			my $label = (grep {(lc $column eq lc $_) || (lc "$column\_0" eq lc $_) ||
					       (lc $column eq lc "$_\_0")} @labels)[0];

			if ($label) {
			    $column = $label; # the column name from the text file is used to match the actual data
			}
			else {
			    ($outfile eq '>') || 
				warn ("Warning: could not find column '$column' in table '$table', file \"$file\" at line $.!\n");
			}
			$values{$column} = ''; # initialise (will stay empty for missing columns and if have '---')

		    }

		}
		else {

		    # blank line signifies end of relevant table so jump out of loop
		    last COLUMNS;

		}
		@output = ();
	    }

	      # hit end of file, either successfully or not

	      unless ($tablename) {
		  # For solid phases, or any phase in the trace elements output, a whole table could be missing
		  ($outfile eq '>') || warn ("Warning: could not find table '$table' in file \"$file\" at line $.!\n");
		  my $i = @columns; # number of columns
		  if ($header eq 'matlab') { 
		      $lines[1] .= (join "_$table$delim", @columns)."_$table$delim";
		  }
		  else {
		      $lines[0] .= "\"$table\"";
		      $lines[1] .= (join $delim, @columns).$delim;
		  }
		  $lines[0] .= $delim x $i;
	      }

	      while($lineno < @pt1) { # length of @pt1 array is zero if no variables are set
		  my $i = @columns; # nuber of columns
		  $lines[$lineno++] .= $delim x $i; # else pad the rest of the printed table
	      }
	      ($outfile eq '>') || warn("Processed columns: @columns\n\n");

	      close FILE;

	      $lineno = 2;
	      $tablename = ''; # reset for next go
	      @output = ();
	      %values = ();
	      @columns = ();
	      @labels = ();
	      
	  }

      }

	# With incorrect line endings the whole file may be one long line
	if ($. < 3) {
	    ($outfile eq '>') || warn ("Please run file_format.command on \"$colfile\" and try again!\n");
	    last FOREACH;
	}
	close COLUMNLIST;

    }

}

# Remove extra header line in matlab version
shift @lines if ($header eq 'matlab');

# Output with whatever the correct line endings for the system are
if (@ARGV) {

    print map {$_ .= "\n"} @lines;
    ($outfile eq '>') || warn "Wrote output.\n";

}
else { 

    $outfile =~ s/^(\"|\')//; # Trim any leftover quotes from start and end of file name
    $outfile =~ s/(\"|\')$//;
    # Fix any apostrophes in the name (if had single quotes or no quotes originally)
    $outfile =~ s/\'\\\'\'/\'/g || $outfile =~ s/\\\'/\'/g;
    $outfile =~ s/\\//g unless $windows; # Remove escaping from any other characters (e.g. '(' and ')')

    if (open(OUTPUT, ">$outfile")) {
	print OUTPUT map {$_ .= "\n"} @lines;
	close (OUTPUT);
	print ("Wrote output.\n");
    }
    else {
	warn "Error: could not open file \"$outfile\" : $!\n";
    }
    print ("Press return to finish.\n");
    $_ = <STDIN>;

}
