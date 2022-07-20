#!/usr/local/bin/perl
use Getopt::Long;
use strict;
use warnings;
my ($taxonbase, $taxontxt, $projectname, $type, $basehref, $tt);
&GetOptions("d=s"       => \$taxonbase,
	    "i=s" =>\$taxontxt,
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

unless (-d $taxonbase){
    mkdir($taxonbase, 0755) || die $!;
}

my $bin_dir = "/export/home/xdong/bin";
my $jsonfile = "$taxonbase/tree.json";
my $cmd = "$^X $bin_dir/taxons2jsontree -i $taxontxt -d 3 > $jsonfile";
system($cmd);

my $taxontree = "$taxonbase/taxonTree.html";
print_taxonTree($projectname, $type, $taxontree, $basehref);

my $taxonjs = "$taxonbase/taxons.js";
print_taxonjs($projectname, $type, $taxonjs, $basehref);
chmod(0755, $jsonfile) or die "could not chang $jsonfile permission, $!\n";
chmod(0755, $taxonjs) or die "could not chang $taxonjs permission, $!\n";



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
    <h3>$projectname Interactive Taxonomic Tree built from $tt </h3> 
  
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

