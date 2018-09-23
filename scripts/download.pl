#!/usr/bin/perl

	
	my $accn = shift;

	next unless ($accn =~ m/ERR/);

	print "$accn\n";

	my $six = substr($accn, 0, 6);

	my $asp = `which ascp`;
	chomp($asp);
	$asp =~ s!bin/ascp!etc/asperaweb_id_dsa.openssh!;


	for($i=0;$i<=9;$i++) {

		print "$accn\n";

		my $vol = "00$i";

		my $url = "/vol1/fastq/$six/$vol/$accn/${accn}_1.fastq.gz";


		print "ascp -QT -l 300m -P33001 -i $asp era-fasp\@fasp.sra.ebi.ac.uk:$url ./\n";
		system("ascp -QT -l 300m -P33001 -i $asp era-fasp\@fasp.sra.ebi.ac.uk:$url ./");


		my $url = "/vol1/fastq/$six/$vol/$accn/${accn}_2.fastq.gz";

		system("ascp -QT -l 300m -P33001 -i $asp era-fasp\@fasp.sra.ebi.ac.uk:$url ./");

		last if (-f "${accn}_1.fastq.gz" && -f "${accn}_2.fastq.gz");
	}
