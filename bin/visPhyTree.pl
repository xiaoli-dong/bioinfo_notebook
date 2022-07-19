#!/usr/local/bin/perl

BEGIN{
  defined $ENV{MAGPIEHOME} or 
    die "MAGPIEHOME is not a defined environment variable\n";
  push @INC, "$ENV{MAGPIEHOME}/lib";
}

use Taxonomy;
use Getopt::Long;
use strict;



my ($taxonbase, $stats_bin, $projectname, $type, $basehref, $tt);
&GetOptions("d=s"       => \$taxonbase,
            "i=s"       => \$stats_bin,
	    "p=s" =>\$projectname,
	    "s=s" =>\$type,
	    "t=s" =>\$tt,
	    "h=s" =>\$basehref
    );

($projectname) ||
    die "usage: $0 OPTIONS

where options are:  -d  <base dir>
                    -i  <input file which is the output from SOrt-items>
-p <project name>
-s <454|sanger|illumina>
-t <free text display at the top of the page>
-h <base href>
\n";

my $taxontxt = "$taxonbase/taxons.txt";

unless (-d $taxonbase){
    mkdir($taxonbase, 0755) || die $!;
}
open (OUT, ">$taxontxt")  or die "Can't open $taxontxt for writting: $1\n";

my %taxons = ();
my %stats = ();
my %b_phylum = ();
my %a_phylum = ();
my %superkingdom = ();

open(STATS, "$stats_bin") or die "Could not open $stats_bin to read, $!\n";

while (<STATS>){
    
    #2       Bacteria        5287
    if(/^#/){

	next;
    }
    elsif(/^(\d+?)\t+(.*?)\t(\d+)/){
	my $taxid = $1;
	if (defined $taxid) {
	    my $lineage = taxid2taxon($taxid, $3);
	    if($lineage eq ""){
		$taxons{"$1\_$2"} += $3;
	    }
	    else{
		$taxons{taxid2taxon($taxid, $3)} += $3;
	    }
	}
    }
    elsif(/^(NHBin|UnAss|TAss|TReads)/){
	if(/^(\w+?)\s+.*?(\d+)$/){
	    $stats{$1} = $2;
	}
    }
}

foreach (sort {$taxons{$b} <=> $taxons{$a}} keys %taxons){
    print OUT "$_\t", $taxons{$_}, "\n";
}
close(OUT);



my $jsonfile = "$taxonbase/tree.json";

my $cmd = "$^X /export/home/xdong/bin/taxons2jsontree -i $taxontxt >$jsonfile";
system($cmd);
my $taxonpiechart = "$taxonbase/index.html";
print_taxonPiechart($projectname, $type, $taxonpiechart, $basehref);

my $taxontree = "$taxonbase/taxonTree.html";
print_taxonTree($projectname, $type, $taxontree, $basehref);

my $taxonjs = "$taxonbase/taxons.js";
print_taxonjs($projectname, $type, $taxonjs, $basehref);

my $highchartJS = "$taxonbase/highchartJS.js";
print_highchart_javaScripts($highchartJS, \%stats, \%superkingdom, \%b_phylum, \%a_phylum);

chmod(0755, $taxontxt) or die "could not chang $taxontxt permission, $!\n";
chmod(0755, $jsonfile) or die "could not chang $jsonfile permission, $!\n";
chmod(0755, $taxonpiechart) or die "could not chang $taxonpiechart permission, $!\n";
chmod(0755, $taxontree) or die "could not chang $taxonpiechart permission, $!\n";
chmod(0755, $taxonjs) or die "could not chang $taxonjs permission, $!\n";


sub taxid2taxon{

    my ($taxida, $count) = @_;
    my $tax = new Taxonomy();
    my @taxids = ();
    @taxids = $tax->lineage($taxida);
    my $taxon = "";
    
    if(@taxids == 0){
	#print STDERR "there is no linearge for taxida=$taxida\n";
	return "";
    }
    
    foreach(@taxids){
	my $rank = $tax->rank($_);

	my $sysname = $tax->taxid2sysname($_);

	if($rank eq "superkingdom"){
	    $superkingdom{$sysname} += $count; 
	}elsif($rank eq "phylum" && $taxon =~ /Bacteria/){
	    $b_phylum{$sysname} += $count; 
	}elsif($rank eq "phylum" && $taxon =~ /Archaea/){
	    $a_phylum{$sysname} += $count; 
	}
	$taxon .= "$_\_$sysname;";
    }
    #print STDERR "$taxon\n";
    return $taxon;
}
sub gi2taxon{
    
    my ($gi) = @_;
    
    my $tax = new Taxonomy();
    my $taxida = $tax->gid2taxid($gi);
    my @taxids = ();
    @taxids = $tax->lineage($taxida);
    my $taxon = "";
    foreach(@taxids){
	my $rank = $tax->rank($_);
	my $sysname = $tax->taxid2sysname($_);
	$taxon .= "$sysname;";
    }
    return $taxon;
}
sub print_taxonPiechart{

    my ($projectname, $type, $taxonpiechart, $basehref) = @_;
    open (OUT, ">$taxonpiechart") or die "Could not open $taxonpiechart file to wrtie, $!\n";
    
print OUT <<END;
<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>Taxonomy Tree</title>
       <link rel="stylesheet" type="text/css" href="/HMP/extjs/resources/css/ext-all.css">
       
	<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>
	<script type="text/javascript" src="/HMP/Highcharts/js/highcharts.js"></script>
	<script type="text/javascript" src="/HMP/Highcharts/js/modules/exporting.js"></script>
	<script type="text/javascript" src="/HMP/metagenomes/data/$projectname/$type/$taxonbase/highchartJS.js"></script>


        <style type="text/css">

	body {
	    font-size:100%;
	    margin :20;
    }
    h1 {font-size:2.0em;}
    h2 {font-size:1.875em;}
    h3 {font-size:1.15em;}
    p {font-size:0.875em;}
          
    
</style>
    </head>
    <body>
     
    <h1>Taxonomic Classification Summaries: $tt  </h1><br><br>
    <table>
    <tr><td id="stats"></td><td id="superkingdom" ></td></tr>
    <tr><td id="bphylum"></td><td id="aphylum"></td></tr></table>
    
    </body>
</html>
END
    
close(OUT);

}


sub print_taxonTree{

    my ($projectname, $type, $taxontree, $basehref) = @_;
    open (OUT, ">$taxontree") or die "Could not open $taxontree file to wrtie, $!\n";
    
print OUT <<END;
<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
        <title>Taxonomy Tree</title>
        <link rel="stylesheet" type="text/css" href="/HMP/extjs/resources/css/ext-all.css">
	<link rel="stylesheet" type="text/css" href="/HMP/extjs/resources/css/mask.css">
	<script type="text/javascript" src="/HMP/extjs/bootstrap.js"></script>
        <script type="text/javascript" src="/HMP/$basehref/$taxonbase/taxons.js"></script>

        <style type="text/css">
	.task {
                background-image: url(/HMP/extjs/examples/shared/icons/fam/cog.gif) !important;
    }
    .task-folder {
	background-image: url(/HMP/extjs/examples/shared/icons/fam/folder_go.gif) !important;
    }
body {font-size:100%;}
    h1 {font-size:2.0em;}
    h2 {font-size:1.875em;}
    h3 {font-size:1.15em;}
p {font-size:0.875em;}
          
    #tree-example{
    
  background: #FFFFFF;
    //border: 2px solid #ccc;
	//margin-top: 50px;
    margin-left: 10px;
  padding: 10px 15px;
  color: #000;
    
}
#navigation {
 position: absolute;
 top: 0;
 left: 0;
 width: 10em;
}

#content {
margin-top: 10px;
margin-left: 10px;
padding: 10px 15px;
  color: #000;
}


</style>
    </head>
    <body>
     
    <div id="content">
    <h2>Interactive Taxonomic Tree: $tt </h1> 
  
    </div>
    <div id="loading-mask"></div>
    <div id="loading">
    <div class="loading-indicator">
    Loading...
    </div>
    </div>
    <div id="tree-example"></div>
    </body>
</html>
END
    
close(OUT);

}


sub print_taxonjs{

    my ($projectname, $type, $taxonjs, $basehref) = @_;
    open (OUT, ">$taxonjs") or die "Could not open $taxonjs file to wrtie, $!\n";

print OUT <<END;
/*

This file is part of Ext JS 4

Copyright (c) 2011 Sencha Inc

Contact:  http://www.sencha.com/contact

GNU General Public License Usage
This file may be used under the terms of the GNU General Public License version 3.0 as published by the Free Software Foundation and appearing in the file LICENSE included in the packaging of this file.  Please review the following information to ensure the GNU General Public License version 3.0 requirements will be met: http://www.gnu.org/copyleft/gpl.html.

If you are unsure which license is appropriate for your use, please contact the sales department at http://www.sencha.com/contact.

*/
Ext.require([
    'Ext.data.*',
    'Ext.grid.*',
    'Ext.tree.*'
]);

Ext.onReady(function() {

    //loading mask
   
    setTimeout(function(){
   Ext.get('loading').remove();
  Ext.get('loading-mask').fadeOut({remove:true});
}, 350);

    //we want to setup a model and store instead of using dataUrl
    Ext.define('Task', {
        extend: 'Ext.data.Model',
        fields: [
            {name: 'task',     type: 'string'},
            {name: 'duration', type: 'int'},
	    {name: 'percent', type: 'float'}
        ]
    });

    var store = Ext.create('Ext.data.TreeStore', {
        model: 'Task',
        proxy: {
            type: 'ajax',
            //the store will get the content from the .json file
            //url: 'taxontree.json'
            url: 'tree.json'
            //url: 'realtree.json'
        },
        folderSort: true
    });

    //Ext.ux.tree.TreeGrid is no longer a Ux. You can simply use a tree.TreePanel
    var tree = Ext.create('Ext.tree.Panel', {
        title: 'Taxonomy Tree Visualization',
        width: 600,
        height: 800,
        //renderTo: Ext.getBody(),
	  renderTo: 'tree-example',
        collapsible: true,
        useArrows: true,
        rootVisible: false,
        store: store,
        multiSelect: true,
        //singleExpand: true,
	  viewConfig: {
            stripeRows: true
        },
        //the 'columns' property is now 'headers'
        columns: [{
            xtype: 'treecolumn', //this is so we know which column will show the tree
            text: 'Taxonomy Terms',
            flex: 2,
            sortable: true,
            dataIndex: 'task'
        },{
            //we must use the templateheader component so we can use a custom tpl
            xtype: 'gridcolumn',
            text: 'Sequence Counts',
            flex: 1,
            sortable: true,
            dataIndex: 'duration',
            align: 'center'
            //add in the custom tpl for the rows

        },{
            //we must use the templateheader component so we can use a custom tpl
            xtype: 'gridcolumn',
            text: 'Count Percentage',
            flex: 1,
            sortable: true,
            dataIndex: 'percent',
            align: 'center'
            //add in the custom tpl for the rows

        }]
    });
});
END
    
    close(OUT);

 
    
}


sub print_highchart_javaScripts{
    
    my($highchartJS, $stats, $superkingdom, $b_phylum, $a_phylum) = @_;
    
    my $s = "";
    foreach (sort {$superkingdom->{$a} <=> $superkingdom->{$b}} keys %$superkingdom){
	 #['TAss', $stats->{"TAss"}],
	$s .= "[\'$_\', " . $superkingdom->{$_} . "],\n";
	
    }
    my $bp = "";
    foreach (sort {$b_phylum->{$a} <=> $b_phylum->{$b}} keys %$b_phylum){
	 #['TAss', $stats->{"TAss"}],
	$bp .= "[\'$_\', " . $b_phylum->{$_} . "],\n";
    }
    my $ap = "";
    foreach ( sort {$a_phylum->{$a} <=> $a_phylum->{$b}}  keys %$a_phylum){
	 #['TAss', $stats->{"TAss"}],
	$ap .= "[\'$_\', " . $a_phylum->{$_} . "],\n";
    }
    
        
    open (HIGHCHART, ">$highchartJS") or die "Could not open $highchartJS  to wrtie, $!\n";

    print HIGHCHART <<END;
    
    var charts = [];
    
    \$(document).ready(function() {
	var getChartConfig = function(renderId, titletxt  , data, dist) {
	    var config = {};
	    config.chart = {
	      renderTo: renderId,
	      defaultSeriesType: 'pie',
	      //borderColor: '#EBBA95',
	      borderWidth: 1,
	      //width: 250,
	      margin: [30, 100, 30, 30]	  
	    };
	    config.credits = {
	      enabled: false
	    };
	    config.title = {
		floating: true,
		verticalAlign: 'bottom',
		y: -5,
		style: {
		    //color: '#FF00FF',
		    fontWeight: 'bold'
		},
		text: titletxt  
	};
	    
	
	    config.tooltip = {
			formatter: function() {
				return '<b>'+ this.point.name +'</b>: '+ this.y;
			}
		},

	    config.plotOptions = {
	      pie: {
		dataLabels: {
		  distance: dist,
		  //color: 'white',
		  formatter: function() {
		      return '<b>'+ this.point.name +'</b>: '+ Highcharts.numberFormat(this.percentage,2) + '%';
  }
  }
  }
  };
  config.navigation = {
        buttonOptions: {
  align: 'right'
 
  }
    };
		      
		      config.legend = { enabled: true };
		      
		      config.series = data;
		      
		      return config;
		  };
  
		    var data = [{
		      data: [
			  ['TAss', $stats->{"TAss"}],
			  ['NHBin', $stats->{"NHBin"}],
			  ['UnAss', $stats->{"UnAss"}],
			] }];
			  
			  var data_s =  [{
		      data: [
			  $s
			      ] }];
			  var data_bp =  [{
		      data: [
			  $bp
			      ] }];
			   var data_ap =  [{
		      data: [
			  $ap
			       ] }];
			  
		       var title_stats = "Taxonomic Classification Stats";
			  var title_kingdom = "Superkingdom Distribution";
			  var title_bphylum = "Bacteria Phylum Distribution";
			  var title_aphylum = "Archaea Phylum Distribution";
			  var dis = 10;
		    //now, creating a new chart is easy!
				charts.push(new Highcharts.Chart(
							 getChartConfig("stats", title_stats, data, dis)
				    ));

		    charts.push(new Highcharts.Chart(
						     getChartConfig("superkingdom", title_kingdom, data_s, dis)
				    ));

		    charts.push(new Highcharts.Chart(
						     getChartConfig("bphylum", title_bphylum, data_bp, dis)
				    ));
	    charts.push(new Highcharts.Chart(
					     getChartConfig("aphylum", title_aphylum, data_ap, dis)
				    ));
		    
		  });
END
		       close(HIGHCHART);
		       
	}
		  
