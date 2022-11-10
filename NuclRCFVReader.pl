use 5.014;
use warnings;
use FAST::Bio::SeqIO;

print "
RCFV Reader (Nucleotide Version), James F. Fleming & Torsten H Struck, 2022.
Welcome to RCFV Reader, Nucleotide Version. RCFV Reader accepts nucleotide FASTA files as input, and then outputs 6 files.
To run RCFV Reader, input a command as follows:

perl NuclRCFVReader.pl <filename> <prefix for output files>

the 6 output files are:
- RCFV.txt - this file contains the total RCFV and nRCFV of the entire dataset.
- csRCFV.txt - this file contains a list of the character-specific RCFVs of your input dataset
- ncsRCFV.txt - as per csRCFV.txt, but with ncsRCFV values
- tRCFV.txt - this file contains a list of the taxon-specific RCFVs of your input dataset
- ntRCFV.txt - as per tRCFV.txt, but with ntRCFV values
- Frequencies.txt - this file contains the per taxon character frequencies, and the mean per taxon character frequencies across the dataset (as used by RCFV, NOT the total dataset character frequencies)

I hope you enjoy this fast new way to use normalised Relative Frequency Composition Values! If you have any questions or queries, or notice any bugs, contact me at:
j.f.fleming\@nhm.uio.no
";

my $fasta  = FAST::Bio::SeqIO->new(-file => $ARGV[0], -format => 'Fasta');
my $seqnum=0;
my $RCFV_File = "$ARGV[1].RCFV.txt";
my $csRCFV_File = "$ARGV[1].csRCFV.txt";
my $tRCFV_File = "$ARGV[1].tRCFV.txt";
my $ncsRCFV_File = "$ARGV[1].ncsRCFV.txt";
my $ntRCFV_File = "$ARGV[1].ntRCFV.txt";
my $outfile = "$ARGV[1].Frequencies.txt";
open(OUT, '>', $outfile) or die $!;
open(RCFV_OUT, '>', $RCFV_File) or die $!;
open(CSRCFV_OUT, '>', $csRCFV_File) or die $!;
open(TRCFV_OUT, '>', $tRCFV_File) or die $!;
open(NCSRCFV_OUT,  '>', $ncsRCFV_File) or die $!;
open(NTRCFV_OUT, '>', $ntRCFV_File) or die $!;

my %A_Lens;
my %C_Lens;
my %G_Lens;
my %T_Lens;
my %A_Freqs;
my %C_Freqs;
my %G_Freqs;
my %T_Freqs;
my @ID_List;
my @Length;

print OUT "NAME\tFreq(A)\tFreq(C)\tFreq(G)\tFreq(T)\n";
while ( my $seq = $fasta->next_seq() ) {
    my $stats;
    $stats->{len} = length($seq->seq);
    push (@Length, length($seq->seq));
    $stats->{$_}++ for split //, $seq->seq;
 #   say ++$seqnum, " @$stats{qw(len A C G T)}";
    my $seqTotal = @$stats{qw(A)} + @$stats{qw(C)} + @$stats{qw(T)} + @$stats{qw(G)};
    my $seqA_Freq = @$stats{qw(A)}/$seqTotal;
    my $seqC_Freq = @$stats{qw(C)}/$seqTotal;
    my $seqG_Freq = @$stats{qw(G)}/$seqTotal;
    my $seqT_Freq = @$stats{qw(T)}/$seqTotal;
	print OUT $seq->id, "\t", $seqA_Freq, "\t", $seqC_Freq, "\t", $seqG_Freq, "\t", $seqT_Freq, "\n";
	push (@ID_List, ($seq->id));
    $A_Lens{$seq->id} = @$stats{qw(A)};
	$C_Lens{$seq->id} = @$stats{qw(C)};
	$G_Lens{$seq->id} = @$stats{qw(G)};
	$T_Lens{$seq->id} = @$stats{qw(T)};
    $A_Freqs{$seq->id} = $seqA_Freq;
    $C_Freqs{$seq->id} = $seqC_Freq;
    $G_Freqs{$seq->id} = $seqG_Freq;
    $T_Freqs{$seq->id} = $seqT_Freq;
}

my %check;
@check{@Length} = (1) x @Length;
if (keys %check == 1){
	print "Fasta File Aligned\n";
	} 
	else{
	print "Are you sure this file is aligned? Sequences seem to have differing lengths.\n" and die $!;}


my $A_total = eval join '+', values %A_Freqs;
my $C_total = eval join '+', values %C_Freqs;
my $G_total = eval join '+', values %G_Freqs;
my $T_total = eval join '+', values %T_Freqs;
my $Total = $A_total + $C_total + $G_total + $T_total;

#print $Total, "\t", $A_total, "\t", $C_total, "\t", $G_total, "\t", $T_total, "\n";

my $MeanA_Freq = ($A_total/$Total);
my $MeanC_Freq = ($C_total/$Total);
my $MeanG_Freq = ($G_total/$Total);
my $MeanT_Freq = ($T_total/$Total);

print OUT "Mean_Freq_Across_Taxa\t", $MeanA_Freq, "\t", $MeanC_Freq, "\t", $MeanG_Freq, "\t", $MeanT_Freq, "\n";

my %A_RCFV_List;
my %C_RCFV_List;
my %G_RCFV_List;
my %T_RCFV_List;

foreach my $Akey (keys %A_Freqs){
	my $Mu =  abs($A_Freqs{$Akey}-$MeanA_Freq);
	$A_RCFV_List{$Akey} = $Mu;
#	print $Mu, "\n";
	}
my $A_Size = keys %A_RCFV_List;
my $A_RCFV_Total = eval join '+', values %A_RCFV_List;
my $A_RCFV = $A_RCFV_Total/$A_Size;

foreach my $Ckey (keys %C_Freqs){
	my $Mu =  abs($C_Freqs{$Ckey}-$MeanC_Freq);
	$C_RCFV_List{$Ckey} = $Mu;
#	print $Mu, "\n";
	}
my $C_RCFV_Total = eval join '+', values %C_RCFV_List;
my $C_RCFV = $C_RCFV_Total/$A_Size;

foreach my $Gkey (keys %G_Freqs){
	my $Mu =  abs($G_Freqs{$Gkey}-$MeanG_Freq);
	$G_RCFV_List{$Gkey} = $Mu;
#	print $Mu, "\n";
	}
my $G_RCFV_Total = eval join '+', values %G_RCFV_List;
my $G_RCFV = $G_RCFV_Total/$A_Size;

foreach my $Tkey (keys %T_Freqs){
	my $Mu =  abs($T_Freqs{$Tkey}-$MeanT_Freq);
	$T_RCFV_List{$Tkey} = $Mu;
#	print $Mu, "\n";
	}
my $T_RCFV_Total = eval join '+', values %T_RCFV_List;
my $T_RCFV = $T_RCFV_Total/$A_Size;

print TRCFV_OUT "\ntRCFV values:\n";
print NTRCFV_OUT "\nntRCFV values:\n";

my $i= 0;
my $tRCFV;
foreach(@ID_List){
	my $unlock=$_;
	$tRCFV = ($A_RCFV_List{$unlock}+$C_RCFV_List{$unlock}+$G_RCFV_List{$unlock}+$T_RCFV_List{$unlock})/$A_Size;
	print TRCFV_OUT $_,"\t",$tRCFV,"\n";
	print NTRCFV_OUT $_,"\t",(0.25*$tRCFV*$A_Size*sqrt($Length[0]))/100,"\n";
	$i++;
	}

my $RCFV = $A_RCFV + $C_RCFV + $G_RCFV + $T_RCFV;
my $nRCFV = $RCFV/(($Length[0]**-0.5)*400);

print CSRCFV_OUT "\ncharacter RCFV values:\ncsRCFV(A)\t", $A_RCFV, "\ncsRCFV(C)\t", $C_RCFV, "\ncsRCFV(G)\t", $G_RCFV, "\ncsRCFV(T)\t", $T_RCFV, "\n";

print NCSRCFV_OUT "\ncharacter RCFV values:\ncsRCFV(A)\t", $A_RCFV/(($Length[0]**-0.5)*100), "\ncsRCFV(C)\t", $C_RCFV/(($Length[0]**-0.5)*100), "\ncsRCFV(G)\t", $G_RCFV/(($Length[0]**-0.5)*100), "\ncsRCFV(T)\t", $T_RCFV/(($Length[0]**-0.5)*100), "\n";

print RCFV_OUT "\ntotal RCFV\nRCFV\t", $RCFV, "\nnRCFV\t", $nRCFV, "\n";


&end;

sub end{
		
		TIMER:
		# set timer
		my ( $user, $system, $cuser, $csystem ) = times ;
		
		print "\n\n\t-------------------------------------------\n\tRCFV CALCULATION COMPLETE Ta'ra\n\t-------------------------------------------\n\t";
		
		print <<TIME;
		
		***  time used: $user sec  ***
		
TIME
		

		
		
		
}
