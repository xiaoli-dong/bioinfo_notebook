use Net::FTP;

#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/995/GCF_000195995.1_ASM19599v1/
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/$str/*/*_protein.faa.gz

my ( $host, $dir ) =
  qw(ftp.ncbi.nih.gov /genomes/all/GCF/000/195/995/GCF_000195995.1_ASM19599v1);

$ftp = Net::FTP->new( $host, Debug => 0 )
  or die "Cannot connect to some.host.name: $@";

$ftp->login( "anonymous", '-anonymous@' ) or die "Cannot login ", $ftp->message;
$ftp->cwd($dir) or die "Cannot change working directory ", $ftp->message;

my @array = grep /_protein\.faa\.gz/, $ftp->ls("$dir/*")
  or die "Cannot ls working directory\n", $ftp->message;

print join( "\n", @array );
print "Fetching remote file $file...\n";

#$ftp->get($file, "$localdir/$file");
use File::Basename;
my $path = $array[0];
print "path=$path\n";
my $filename = basename($path);

print "filename=$filename\n";
$ftp->get( $path, "./$filename" ) or die "get $path failed ", $ftp->message;

$ftp->quit;
