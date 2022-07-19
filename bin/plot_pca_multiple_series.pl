#!/usr/local/bin/perl

$#ARGV == -1

  and die
"Usage: $0 PCA clustering output from fast unifrac, PCA plot template, sample category, title\n";

my $pca_file    = $ARGV[0];
my $template    = $ARGV[1];
my $sample_list = $ARGV[2];
my $title       = $ARGV[3];

my %color_map = ();

open( SAMPLE, "<$sample_list" ) or die $!;

#235     019P3P_454AMP   SuncorTailings Pond 5 (45 ft)   LG_SUN_TP_TS_010_014_071009     TP-E

my %sample_des       = ();
my %name_mapping     = ();
my %category_mapping = ();
while ( my $entry = <SAMPLE> ) {

    if ( $entry =~ /Won't be sequenced/ ) {

        #do nothing
    }
    elsif ( $entry =~ /No Read!!!/ ) {

    }
    elsif ( my ( $GQ_name, $des, $unique_id, $category ) =
        $entry =~ /^\d+\s+(\S+)\s+(\S+.*?)\t+(\S+)\t(\S+)/ )
    {
        print STDERR "category=$category\n";
        $GQ_name =~ s/-/_/g;

        $des =~ s/'/\\'/g;
        $sample_des{"$GQ_name\_$unique_id"}       = $des;
        $name_mapping{"$GQ_name\_$unique_id"}     = $unique_id;
        $category_mapping{"$GQ_name\_$unique_id"} = $category;

        #print STDERR "$gqname\_$id => $des\n";
    }
}

close(SAMPLE);

foreach ( sort keys %sample_des ) {

    #print STDERR "$_ => $sample_des{$_}\n";

}
open( PCA, "<$pca_file" ) or die $!;

my $p1_p2      = "";
my $p1_p3      = "";
my $p1_p4      = "";
my $p2_p3      = "";
my $p2_p4      = "";
my $p3_p4      = "";
my $var1       = "";
my $var2       = "";
my $var3       = "";
my $var4       = "";
my $p1_b       = 0;
my $p2_b       = 0;
my $p3_b       = 0;
my $p4_b       = 0;
my %dataSeries = ();

while ( my $entry = <PCA> ) {
    if ( $entry =~ /var explained \(\%\)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
        $var1 = $1;
        $var2 = $2;
        $var3 = $3;
        $var4 = $4;

        #print STDERR "$var1";
    }
    elsif ( my ( $sampleid, $p1, $p2, $p3, $p4 ) =
        $entry =~ /^Eigenvectors\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ )
    {

        if ( $p1_b < abs($p1) ) {

            $p1_b = abs($p1);
        }
        if ( $p2_b < abs($p2) ) {

            $p2_b = abs($p2);
        }
        if ( $p3_b < abs($p3) ) {

            $p3_b = abs($p3);
        }
        if ( $p4_b < abs($p4) ) {

            $p4_b = abs($p4);
        }

        if ( exists $dataSeries{"12"}->{ $category_mapping{$sampleid} } ) {

            $dataSeries{"12"}->{ $category_mapping{$sampleid} } .=
                ",\n{\nx: $p1,\ny: $p2,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }
        else {
            $dataSeries{"12"}->{ $category_mapping{$sampleid} } =
                "{\nx: $p1,\ny: $p2,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }

        if ( exists $dataSeries{"13"}->{ $category_mapping{$sampleid} } ) {

            $dataSeries{"13"}->{ $category_mapping{$sampleid} } .=
                ",\n{\nx: $p1,\ny: $p3,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }
        else {
            $dataSeries{"13"}->{ $category_mapping{$sampleid} } =
                "{\nx: $p1,\ny: $p3,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }

        if ( exists $dataSeries{"14"}->{ $category_mapping{$sampleid} } ) {

            $dataSeries{"14"}->{ $category_mapping{$sampleid} } .=
                ",\n{\nx: $p1,\ny: $p4,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }
        else {
            $dataSeries{"14"}->{ $category_mapping{$sampleid} } =
                "{\nx: $p1,\ny: $p4,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }

        if ( exists $dataSeries{"23"}->{ $category_mapping{$sampleid} } ) {

            $dataSeries{"23"}->{ $category_mapping{$sampleid} } .=
                ",\n{\nx: $p2,\ny: $p3,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }
        else {
            $dataSeries{"23"}->{ $category_mapping{$sampleid} } =
                "{\nx: $p2,\ny: $p3,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }

        if ( exists $dataSeries{"24"}->{ $category_mapping{$sampleid} } ) {

            $dataSeries{"24"}->{ $category_mapping{$sampleid} } .=
                ",\n{\nx: $p2,\ny: $p4,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }
        else {
            $dataSeries{"24"}->{ $category_mapping{$sampleid} } =
                "{\nx: $p2,\ny: $p4,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }

        if ( exists $dataSeries{"34"}->{ $category_mapping{$sampleid} } ) {

            $dataSeries{"34"}->{ $category_mapping{$sampleid} } .=
                ",\n{\nx: $p3,\ny: $p4,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }
        else {
            $dataSeries{"34"}->{ $category_mapping{$sampleid} } =
                "{\nx: $p3,\ny: $p4,\nsampleName: \'"
              . $name_mapping{$sampleid}
              . "\',\nalarms_tooltip: \'"
              . $sample_des{$sampleid} . "\'\n}";
        }

    }

}
close(PCA);
$p1_b += 0.1;
$p2_b += 0.1;
$p3_b += 0.1;
$p4_b += 0.1;

#print "$p1_p2";

open( TEMPLATE, "<$template" ) or die $!;

my @color = ( "blue", "green", "hotpink", "purple", "orange", "red" );

my $count = 0;

foreach my $cate ( keys %{ $dataSeries{"12"} } ) {

    if ($p1_p2) {
        $p1_p2 .= ",{\nname: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"12"}->{$cate} . "\n]}";

    }
    else {

        $p1_p2 = "{\nname: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"12"}->{$cate} . "\n]}";
    }
    $count++;

}

$count = 0;

foreach my $cate ( keys %{ $dataSeries{"13"} } ) {

    if ($p1_p3) {
        $p1_p3 .= ",{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"13"}->{$cate} . "]}";

    }
    else {

        $p1_p3 = "{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"13"}->{$cate} . "\n]}";
    }
    $count++;
}
$count = 0;

foreach my $cate ( keys %{ $dataSeries{"14"} } ) {

    if ($p1_p4) {
        $p1_p4 .= ",{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"14"}->{$cate} . "\n]}";

    }
    else {

        $p1_p4 = "{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"14"}->{$cate} . "\n]}";
    }
    $count++;
}

$count = 0;
foreach my $cate ( keys %{ $dataSeries{"23"} } ) {

    if ($p2_p3) {
        $p2_p3 .= ",{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n"
          . $dataSeries{"23"}->{$cate} . "\n]}";

    }
    else {

        $p2_p3 = "{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"23"}->{$cate} . "\n]}";
    }
    $count++;
}
$count = 0;

foreach my $cate ( keys %{ $dataSeries{"24"} } ) {

    if ($p2_p4) {
        $p2_p4 .= ",{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n"
          . $dataSeries{"24"}->{$cate} . "\n]}";

    }
    else {

        $p2_p4 = "{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"24"}->{$cate} . "\n]}";
    }
    $count++;
}
$count = 0;

foreach my $cate ( keys %{ $dataSeries{"34"} } ) {

    if ($p3_p4) {
        $p3_p4 .= ",{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"34"}->{$cate} . "\n]}";

    }
    else {

        $p3_p4 = "{name: \'$cate\',\n color: \'$color[$count]\',\ndata:[\n "
          . $dataSeries{"34"}->{$cate} . "\n]}";
    }
    $count++;
}

while ( my $entry = <TEMPLATE> ) {

    print "$entry";
}

plot(
    "P1",  "P2",  "-$p1_b",      "$p1_b", "-$p2_b", "$p2_b",
    $var1, $var2, "container12", $p1_p2
);
plot(
    "P1",  "P3",  "-$p1_b",      "$p1_b", "-$p3_b", "$p3_b",
    $var1, $var3, "container13", $p1_p3
);
plot(
    "P1",  "P4",  "-$p1_b",      "$p1_b", "-$p4_b", "$p4_b",
    $var1, $var4, "container14", $p1_p4
);
plot(
    "P2",  "P3",  "-$p2_b",      "$p2_b", "-$p3_b", "$p3_b",
    $var2, $var3, "container23", $p2_p3
);
plot(
    "P2",  "P4",  "-$p2_b",      "$p2_b", "-$p4_b", "$p4_b",
    $var2, $var4, "container24", $p2_p4
);
plot(
    "P3",  "P4",  "-$p3_b",      "$p3_b", "-$p4_b", "$p4_b",
    $var3, $var4, "container34", $p3_p4
);

print "</head>\n";
print "<body>\n";
print
"<h3 align=\"center\">All $title sample PCA analysis result scatter plotting</h3>\n";
print "<!-- 3. Add the container -->\n";
print
"<div id=\"container12\" style=\"width: 800px; height: 500px; margin: 0 auto\"></div>\n";
print
"<div id=\"container13\" style=\"width: 800px; height: 500px; margin: 0 auto\"></div>\n";
print
"<div id=\"container14\" style=\"width: 800px; height: 500px; margin: 0 auto\"></div>\n";

print
"<div id=\"container23\" style=\"width: 800px; height: 500px; margin: 0 auto\"></div>\n";
print
"<div id=\"container24\" style=\"width: 800px; height: 500px; margin: 0 auto\"></div>\n";
print
"<div id=\"container34\" style=\"width: 800px; height: 500px; margin: 0 auto\"></div>\n";

print "</body>\n";
print "</html>\n";

sub plot {

    my ( $x, $y, $xmin, $xmax, $ymin, $ymax, $var1, $var2, $container, $data )
      = @_;

    print "<script type=\"text/javascript\">\n";
    print "var chart;\n";
    print "\$(document).ready(function() {\n";
    print "chart = new Highcharts.Chart({\n";
    print "chart: {\n";
    print "renderTo: \'$container\',\n";
    print "defaultSeriesType: \'scatter\',\n";
    print "zoomType: \'xy\'\n";
    print "	},\n";
    print "	title: {\n";
    print "text: \'PCoA - $x vs $y\'\n";
    print "},\n";

    #print "legend: {align: 'left', verticalAlign: 'top'},\n";
    print "credits: {enabled: false},\n";
    print "xAxis: {\n";
    print "plotLines: [{ \n";
    print "color: \'#FF0000\',\n";
    print "width: 2,\n";
    print "value: 0,\n";
    print "id: \'plotline-1\'\n";
    print "}],\n";
    print "title: {\n";
    print "enabled: true,\n";
    print "text: \'$x - Percent variation explained $var1 \%\'\n";
    print "},\n";

    print "min: $xmin,\n";
    print "max: $xmax\n";
    print "},\n";
    print "yAxis: {\n";

    print "plotLines: [{ \n";
    print "color: \'#FF0000\',\n";
    print "width: 2,\n";
    print "value: 0,\n";
    print "id: \'plotline-1\'\n";
    print "}],\n";
    print "title: {\n";
    print "text: \'$y - Percent variation explained $var2 \%\'\n";
    print "},\n";

    print "lineWidth: 2,\n";
    print "min: $ymin,\n";
    print "max: $ymax\n";

    print "},\n";
    print "tooltip: {\n";
    print "formatter: function() {\n";
    print
"return '<strong>' + this.point.sampleName + ':</strong><br/>' + this.point.alarms_tooltip;\n";

    #print "return \'\' + this.point.alarms_tooltip;\n";
    print "}\n";
    print "},\n";

    print "plotOptions: {\n";

    print "scatter: {\n";
    print "marker: {\n";
    print "symbol: 'circle',\n";
    print "radius: 5,\n";
    print "states: {\n";
    print "hover: {\n";
    print "enabled: false,\n";
    print "lineColor: \'rgb(100,100,100)\'\n";
    print "}\n";
    print "}\n";
    print "},\n";
    print "states: {\n";
    print "hover: {\n";
    print "marker: {\n";
    print "enabled: false\n";
    print "}\n";
    print "}\n";
    print "}\n";
    print "}\n";

    print "},\n";

    print "series: [ \n";

    #print "name: '',\n";
    #print "color: 'rgba(119, 152, 191, .5)',\n";
    #print "data: [\n";

    print "$data\n";
    print "]\n";

    #print "}\n";
    #print "]\n";
    print "});\n";

    print "});\n";

    print "</script>\n";

}
