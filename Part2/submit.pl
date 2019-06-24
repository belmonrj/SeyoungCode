#!/usr/bin/perl -w

use File::Copy;
use Cwd;

open (FILE, "<list_good_runnum_seg.lst") || die "File Open Error!!\n";
@runlist = <FILE>;
close(FILE);

$maindir = getcwd();
$jobdir = "${maindir}/jobdir_dAu";

mkdir $jobdir;

foreach $runnum (@runlist){

	chomp($runnum);
	
	$wrkdir = "${jobdir}/run${runnum}";
	mkdir $wrkdir;
	chdir $wrkdir;

	open (FILE, ">condor");
	print FILE "Universe      = vanilla\n";
	print FILE "Notification  = Error\n";
	print FILE "Requirements  = CPU_Speed >= 2\n";
	print FILE "Rank          = CPU_Speed\n";
	print FILE "Priority      = +20\n";
	print FILE "request_memory= 1800M\n";
	print FILE "request_cpus  = 1\n";
	print FILE "Initialdir    = ${wrkdir}\n";
	print FILE "Executable    = jobscript\n";
	print FILE "Output        = jobscript.out\n";
	print FILE "Error         = jobscript.err\n";
	print FILE "Log           = jobscript.log\n";
	print FILE "Notify_user   = seyoung922\@gmail.com\n";
	print FILE "Queue";
	close (FILE);

	open (FILE, ">jobscript");
	print FILE "#!/bin/csh -f\n";
	print FILE "setenv HOME /phenix/u/\$LOGNAME\n";
	print FILE "source /etc/csh.login\n";
	print FILE "foreach i (/etc/profile.d/*.csh)\n";
	print FILE "source \$i\n";
	print FILE "end\n";
	print FILE "source \$HOME/.login\n";
	print FILE "source /opt/phenix/bin/phenix_setup.csh\n";

	print FILE "ln -sf ${maindir}/pAu_Cor.C\n";
	print FILE "ln -sf ${maindir}/pAu_Cor.h\n";
	print FILE "ln -sf ${maindir}/def_Const.h\n";
	print FILE "ln -sf ${maindir}/Run.C\n";

	$infile = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/seyoung/taxi/Run16dAu200MuonsMBP107/13874/data/${runnum}.root";
	print FILE "ln -sf ${infile} infile.root\n";

	print FILE "root -l -b -q 'Run.C'\n";
	$outfile = "${jobdir}/outfile_${runnum}.root";
	print FILE "mv outfile.root ${outfile}\n";

	close (FILE);
	chmod 0755, "jobscript";
	system "condor_submit condor";

}
