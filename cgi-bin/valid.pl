#!/usr/bin/perl

############################### Header Information ##############################
require 'cgi.perl';
use CGI;;
$query = new CGI;
&ReadParse;
print &PrintHeader;

################################ Reads Input Data ##############################
$atom = $query->param('atom');
$file = $query->param('file');
$svm_th = $query->param('svm_th');

#################Validation Of Input Sequence Data (file upload) ###################################
if($file ne '' && $atom eq '')
{
    $file=~m/^.*(\\|\/)(.*)/; 
    while(<$file>) 
    {
	$seqfi .= $_;
    }
}
elsif($atom ne '' && $file eq ''){

    $seqfi="$atom";
}

##############ACTUAL PROCESS BEGINS FROM HERE#######################
$infut_file = "/webservers/cgi-bin/erpred";
$ran= int(rand 10000);
$dir = "/webservers/cgidocs/mkumar/temp/Ravindra/ERPred/erpred$ran";
system "mkdir $dir";
system "chmod 777 $dir";
$nam = 'input.'.'fasta';
open(FP1,">$dir/input_meta.fasta");
print FP1 "$seqfi\n";
#print "$seqfi\n";
close FP1;

system "/usr/bin/tr -d '\r' <$dir/input_meta.fasta >$dir/input_fi.fasta"; #Remove meta-character
system "/usr/bin/perl $infut_file/fasta.pl $dir/input_fi.fasta |/usr/bin/head -50 >$dir/input.fasta"; #Convert two line fasta file and select only 25 sequence
system "/bin/grep -c '>' $dir/input.fasta >$dir/total_seq";
system "/bin/grep '>' $dir/input.fasta |/usr/bin/cut -d '|' -f3 |/usr/bin/cut -d ' ' -f1 >$dir/protein_id"; #Grep protein id
$total_seq=`head -1 $dir/total_seq`;chomp($total_seq);
if($total_seq ne 0)
{
    open (FH,"$dir/input.fasta") or die "$!";
    while($line=<FH>)
    {
	chomp($line);
	if($line=~ m/^>/)
	{
	    $next=<FH>;
	    chomp ($next);
	    $n_ter=substr($next,0,25);#N-terminal residues
	    $remain=substr($next,25,-25);#Remaining residues
	    $c_ter=substr($next,-25);#C-terminal residues
	    open(MYFILE,">$dir/sub") or die "$!";
	}
	#print MYFILE"$line\n$n_ter\n$line\n$remain\n$line\n$c_ter\n";
	print MYFILE"$line\n$n_ter\n$line\n$c_ter\n$line\n$remain\n";
	system "/usr/bin/perl $infut_file/aaseqformat.pl $dir/sub +1 >$dir/comp";
	open(FILE,"$dir/comp") or die "$!";
	$c=0;
	open(SAAC,">>$dir/saac") or die "$!";
	print SAAC "+1 ";
	while($file=<FILE>)
	{
	    chomp($file);
	    @array=split(/\+1/,$file);
	    @array1=split(/ /,$array[1]);
	    for($a=1;$a<=$#array1;$a++)
	    {
		$c++;
		#print "$c:$array1[$a] ";
		#open(SAAC,">>$dir/saac") or die "$!";
		print SAAC "$c:$array1[$a] ";
	    }
	}
	print SAAC "\n";
	close SAAC;
    }
    system "/usr/local/bin/svm_classify $dir/saac $infut_file/Models/model_er $dir/svm_score_er >/dev/null";
    system "/usr/bin/paste $dir/protein_id $dir/saac $dir/svm_score_er |/usr/bin/tr '\t' '#' >$dir/final";
    open(FINAL,"$dir/final") or die "$!";
    while($line_fi=<FINAL>)
    {
	chomp($line_fi);
	@svm=split(/\#/,$line_fi);
	if($svm[2] < $svm_th)#Non-ER
	{
	    open(PREDICTION,">>$dir/Prediction") or die "$!";
	    print PREDICTION "$svm[0]\tNon-ER\n";
	    close PREDICTION;
	}
	else
	{
	    open(PREDICTION,">>$dir/Prediction") or die "$!";
	    print PREDICTION "$svm[0]\tER-Protein\n";
	    close PREDICTION;
	}
    }
    close FINAL;
    print  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
    print  "<html><HEAD>\n";
    print  "<TITLE>ERPred::Prediction Result</TITLE>\n";
    print  "<META NAME=\"description\" CONTENT=\"ERPred, University of Delhi South Campus, INDIA\">\n";
    print  "</HEAD><body bgcolor=\"\#FFFFE0\">\n";
    print  "<h2 ALIGN = \"CENTER\"> ERPred Prediction Result</h2>\n";
    print  "<HR ALIGN =\"CENTER\"> </HR>\n";
    print  "<p align=\"center\"><font size=4 color=black><b>The submitted protein/proteins belongs to <font color='red'></p>";
    print "<table border='1' width='400' align='center'><tr><th>Protein ID</th><th>Prediction</th></tr>";
}
if($total_seq == 0)
{
    system "/bin/cat $dir/Pfam_result $dir/svm_result >$dir/Prediction";
    print  "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
    print  "<html><HEAD>\n";
    print  "<TITLE>ERPred::Prediction Result</TITLE>\n";
    print  "<META NAME=\"description\" CONTENT=\"ERPred, University of Delhi South Campus, INDIA\">\n";
    print  "</HEAD><body bgcolor=\"\#FFFFE0\">\n";
    print  "<h2 ALIGN = \"CENTER\"> ERPred Prediction Result</h2>\n";
    print  "<HR ALIGN =\"CENTER\"> </HR>\n";
    print  "<p align=\"center\"><font size=4 color=black><b>Please submit your sequence in fasta format</b></p>";
}
open(PRED,"$dir/Prediction") or die "$!";
while($pre=<PRED>)
{
    chomp($pre);
    @pred=split(/\t/,$pre);
    print "<tr align='center'><td>$pred[0]</td><td>$pred[1]</td></tr>";
}
print "</table>";
print "</font></b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>Thanks for using ERPred Prediction Server</b></font></p>\n";
print  "<p align=\"center\"><font size=3 color=black><b>If you have any problem or suggestions please contact <a href='mailto:manish@south.du.ac.in'>Dr. Manish Kumar</a></b></font>. Please mention your job number in any communication.</p></br>\n";
print  "<p ALIGN=\"CENTER\"><b>Your job number is <font color=\"red\">$ran</b></font></p>\n";
print  "</body>\n";
print  "</html>\n";
system "chmod 000 $dir";
