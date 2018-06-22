package myFootLog;

use strict; use warnings; use Cwd qw(abs_path); use File::Basename qw(dirname fileparse);
my ($N,$B,$BL,$BU,$GN,$YW,$CY,$RD,$PR,$BR,$GR,$WT,$DGR,$LBU,$LGN,$LCY,$LRD,$LPR) = ("\e[0m","\e[1m","\e[0;30m","\e[0;34m","\e[0;32m","\e[1;33m","\e[0;36m","\e[0;31m","\e[0;35m","\e[0;33m","\e[0;37m","\e[0m",,"\e[1;30m","\e[1;34m","\e[1;32m","\e[1;36m","\e[1;31m","\e[1;35m");
use vars qw($N $B $BL $BU $GN $YW $CY $RD $PR $BR $GR $WT $DGR $LBU $LGN $LCY $LRD $LPR);

BEGIN {
	my $homedir = $ENV{"HOME"};
	my $GITFOLDER;
	if (-e "$homedir/.footLoopPath") {
		($GITFOLDER) = `cat $homedir/.footLoopPath`;
	}
	if (not -e "$homedir/.footLoopPath" or (defined $GITFOLDER and not -d $GITFOLDER)) {
		($GITFOLDER) = dirname(dirname abs_path $0) . '/footLoop/';
		print "\n\n----------------------\nInitializing footLoopPath at $homedir/.footLoopPath!\n----------------------\n\n";
		open (my $out, ">", "$homedir/.footLoopPath") or die "Failed to write to $homedir/.footLoopPath: $!\n";
		print $out "$GITFOLDER";
		close $out;
	}
   my $libPath = $GITFOLDER . "lib/";

	if (not defined $libPath or not defined $GITFOLDER or not -d $GITFOLDER) {
		print "\n\nmyFootLog.pm: BEGIN failed trying to run:\n- script=\e[0;32m$0\e[0m\n- libpath=\e[0;33m$libPath\e[0m\n\n";
		exit 1;
	}
   push(@INC, $libPath);
}

my $homedir = $ENV{"HOME"};
my $footLoopScriptDir = dirname(dirname abs_path $0) . "/footLoop";


############
# MAIN LOG #
############

sub MAINLOG {
	my ($RUNSCRIPT, $VALS, $OPTS) = @_;
	
	# Get GIT versions
	my ($GITHASH, $GITCOMMIT, $GITVERSION, $GITDATE, $GITFOLDER, $GITMESSAGES) = getGit(); 

	# Get UUID
	my $RUNUUID = uuid();
	my $RUNDATE = getDate_simple();

	   $RUNSCRIPT = (dirname abs_path $RUNSCRIPT) . "/" . fileparse($RUNSCRIPT);
	my $RUNSCRIPTMD5 = getMD5_simple($RUNSCRIPT);

	my $MAINLOG = "
#gitFolder=$GITFOLDER
#gitCommit=$GITCOMMIT
#gitVersion=$GITVERSION
#gitDate=$GITDATE
##runScript=$RUNSCRIPT
##runScriptMD5=$RUNSCRIPTMD5;
##runUuid=$RUNUUID
##runDate=$RUNDATE
";
	if (defined $OPTS) {
		my $i = 0;
		while ($OPTS =~ /[a-zA-Z0-9]:?/g) {
			my $opt = $&;
			my $val = (not defined $VALS->[$i]) ? $VALS->[$i] : $opt =~ /:$/ ? "-$opt $VALS->[$i]" : "-$opt";
			$val =~ s/://g if defined $val;
			$i++;
			$MAINLOG .= "##runOpts=$val\n" if defined $val;
		}
	}
	return($MAINLOG);
}

sub getMD5_simple {
   my ($file) = @_;
   my ($md5) = `md5sum $file` =~ /^(\w+)\s+/;
   return ($md5);
}

################
# GIT VERSIONS #
################
sub getGit {
	my $gitHash;
	my @gitLog = split(" ", `git --git-dir $footLoopScriptDir/.git log --oneline|head -n 1`);
	my ($gitCommit, $gitVersion, @gitMessages) = @gitLog == 0 ? ("UNK", "UNK", "NA") : @gitLog;
	my ($gitDate) = `git --git-dir $footLoopScriptDir/.git log | head -n 3 |tail -n 1` =~ /^Date:\s+(.+)$/;
	my ($gitMon, $gitNum, $gitYear) = $gitDate =~ /^\w+ (\w+) (\d+) \d+:\d+:\d+ \d\d(\d+) /;
	$gitDate = "$gitYear" . convertMonth($gitMon) . $gitNum;
	my ($gitFolder) = `cat $homedir/.footLoopPath`;
	$gitHash->{gitCommit} = $gitCommit;
	$gitHash->{gitVersion} = $gitVersion;
	$gitHash->{gitDate} = $gitDate;
	$gitHash->{gitFolder} = $gitFolder;
	$gitHash->{gitMessages} = join("", @gitMessages);
	return ($gitHash, $gitCommit, $gitVersion , $gitDate, $gitFolder, join("", @gitMessages));
}

sub getGitCommit {
	my ($gitHash, @arr) = getGit();
	return $gitHash->{gitCommit};
}
sub getGitVersion {
	my ($gitHash, @arr) = getGit();
	return $gitHash->{gitVersion};
}
sub getGitDate {
	my ($gitHash, @arr) = getGit();
	return $gitHash->{gitDate};
}
sub getGitFolder {
	my ($gitHash, @arr) = getGit();
	return $gitHash->{gitFolder};
}
sub getGitMessages {
	my ($gitHash, @arr) = getGit();
	return $gitHash->{gitMessages};
}

##############
# UUID, DATE #
##############

sub convertMonth {
	my ($month) = @_;
	my @month = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	for (my $i = 0; $i < @month; $i++) {
		if ($month[$i] eq $month) {
			$month = $i + 1 < 10 ? 0 . ($i + 1) : $i + 1;
			return $month;
		}
	}
	return -1;
}
sub uuid {
	my ($uuid) = `uuidgen` =~ /^(.+)\n$/;
	return $uuid;
}

sub getDate {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); $year += 1900;
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
   my $date = "$mday $months[$mon] $year $hour:$min:$sec";
   my $timenow = $hour * 3600 + $min * 60 + $sec;
   return($date);
}

sub getDate_simple {
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); $year += 1900;
   my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
   ($year) = $year =~ /^\d\d(\d+)$/;
   ($mon) = $mon + 1 < 10 ? 0 . ($mon+1) : $mon + 1;
   return $year . $mon . $mday;
 #  my $date = "$mday $months[$mon] $year $hour:$min:$sec";
  # my $timenow = $hour * 3600 + $min * 60 + $sec;
  # return($date);
}

sub date {
   my ($add, $color) = @_;
   ($color = $add and undef $add) if defined $add and $add =~ /^(color=|color|col=|col|y|c)/i;
   $add = 0 if not defined $add;
   return ("[$YW" . getDate($add) . "$N]: ") if not defined $color;
   return ("[" . getDate($add) . "]: ") if defined $color;
}

1;
