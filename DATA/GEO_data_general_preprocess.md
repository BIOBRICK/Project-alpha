# GEO Matrix General Preprocess

## Probe id & Gene id

Most of GEO data was collected through Gene Chip tech, thus the row names of the expression will be a little bit difficult to understand. So I have to transform all these probe id into familiar gene name.

Of course, if the gene id is symbol id, that will be good.

```R
setwd("Your_Work_Path")

#If you didn't 	install the GEOquery package, use the script below:
#BiocManager::install('GEOquery') 
library(GEOquery)
library(dplyr)
library(tidyr)

options(stringsAsFactors = F)

GPL6244 <-getGEO('GPL6244',destdir =".")
GPL6244_anno <- Table(GPL6244)

probe2symbol_df <- GPL6244_anno %>% 
  select(ID,gene_assignment) %>% 
  filter(gene_assignment != "---") %>% 
  separate(gene_assignment,c("drop","symbol"),sep="//") %>% 
  select(-drop)

exprSet<-read.table(file = "GSE42872_series_matrix.txt",header = T)
names(exprSet)[1] <- names(probe2symbol_df)[1]
exprSet$probe_id <- as.character(exprSet$probe_id)

exprSet <- exprSet %>% 
  
  inner_join(probe2symbol_df,by="probe_id") %>% #合并探针的信息
  
  select(-probe_id) %>% #去掉多余信息
  
  select(symbol, everything()) %>% #重新排列，
  
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
  
  arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
  
  distinct(symbol,.keep_all = T) %>% # symbol留下第一个
  
  select(-rowMean) %>% #反向选择去除rowMean这一列
  
  tibble::column_to_rownames(colnames(.)[1]) # 把第一列变成行名并删除

write.csv(exprSet, file = 'GSE42872_expression_matrix.csv')

rm(list = ls())

```

**This script isn't my original creation, I get it from internet and thanks to [biotrainee](<http://www.biotrainee.com/>)  very much for this !**

One day, a friend asked me about precessing the expression matrix with some easier way. So, here I am going to intaregate the most common solution for this problem, however it need a time-cost learning precess.

Using perl:
``` Perl
# This is a perl script:
use strict;
use warnings;
print STDERR "gene symbol column number: ";
my $geneSymbolCol=<STDIN>;
chomp($geneSymbolCol);
$geneSymbolCol--;
my $expFile="probeMatrix.txt";
my $gplFile="ann.txt";
my $expFileWF="geneMatrix.txt";
my %hash=();
my @sampleName=();

open(EXP,"$expFile") or die $!;
while(my $exp=<EXP>)
{
	next if ($exp=~/^(\n|\!)/);
	chomp($exp);
	if($.==1)
	{
		my @expArr=split(/\t/,$exp);
		foreach my $singleName (@expArr)
		{
			$singleName=~s/\"//g;
			if($singleName eq 'ID_REF')
			{
				push(@sampleName,$singleName);
			}
			else
			{
				my @singleArr=split(/\_|\./,$singleName);
				push(@sampleName,$singleArr[0]);
			}
		}
	}
	else
	{
		my @expArr=split(/\t/,$exp);
		for(my $i=0;$i<=$#sampleName;$i++)
		{
			$expArr[$i]=~s/\"//g;
			push(@{$hash{$sampleName[$i]}},$expArr[$i]);
		}
	}
}
close(EXP);

my %probeGeneHash=();

open(GPL,"$gplFile") or die $!;
while(my $gpl=<GPL>)
{
	next if($gpl=~/^(\#|ID|\!|\n)/);
	chomp($gpl);
	my @gplArr=split(/\t/,$gpl);
	if((exists $gplArr[$geneSymbolCol]) && ($gplArr[$geneSymbolCol] ne '') && ($gplArr[$geneSymbolCol] !~ /.+\s+.+/))
	{
		$gplArr[$geneSymbolCol]=~s/\//\./g;
		$gplArr[$geneSymbolCol]=~s/\"//g;
		$probeGeneHash{$gplArr[0]}=$gplArr[$geneSymbolCol];
	}
}
close(GPL);

my @probeName=@{$hash{"ID_REF"}};
delete($hash{"ID_REF"});

my %geneListHash=();
my %sampleGeneExpHash=();
foreach my $key (keys %hash)
{
	my %geneAveHash=();
	my %geneCountHash=();
	my %geneSumHash=();
	my @valueArr=@{$hash{$key}};
	for(my $i=0;$i<=$#probeName;$i++)
	{
		if(exists $probeGeneHash{$probeName[$i]})
		{
			my $geneName=$probeGeneHash{$probeName[$i]};
			$geneListHash{$geneName}++;
			$geneCountHash{$geneName}++;
			$geneSumHash{$geneName}+=$valueArr[$i];
		}
	}
	foreach my $countKey (keys %geneCountHash)
	{
		$geneAveHash{$countKey}=$geneSumHash{$countKey}/$geneCountHash{$countKey};
	}
	$sampleGeneExpHash{$key}=\%geneAveHash;
}

open(WF,">$expFileWF") or die $!;
$sampleName[0]="geneNames";
print WF join("\t",@sampleName) . "\n";
foreach my $probeGeneValue (sort(keys %geneListHash))
{
	print WF $probeGeneValue . "\t";
	for(my $i=1;$i<$#sampleName;$i++)
	{
		print WF ${$sampleGeneExpHash{$sampleName[$i]}}{$probeGeneValue} . "\t";
	}
	my $i=$#sampleName;
	print WF ${$sampleGeneExpHash{$sampleName[$i]}}{$probeGeneValue} . "\n";
}
close(WF);
```
In order to use this script, you must have perl. This is not a problem in *nux system in which perl is alread installed. If you a using a windows PC, you need [ActivePerl](https://www.activestate.com/products/activeperl/) or [StrawberryPerl](http://strawberryperl.com), both softwares(You can see these as some kind of software......) are very easy to install.

And notice, there is on row written like this:
```perl
print STDERR "gene symbol column number: ";
```
So you have to enter the number of "Symbol id" in the GPL Platform file downloaded from GEO. And of course an expression matrix is also needed.

Put these two files in the same directory and rename them into "probeMatrix.txt" & "ann.txt", and type perl Your_Script_name.pl in you terminal or PowerShell. After that typing the number 11(The most common number......) and wait, after that a file named "geneMatrix.txt" will be created in the same folder.

But of course this script isn't written by me. And thanks to the unknown one who created this script in the first place!

I have accidently find another method to complete this mission with python! I must cite:

- K. Koizumi, M. Oku, S. Hayashi, A. Inujima, N. Shibahara, L. Chen, Y. Igarashi, K. Tobe, S. Saito, M. Kadowaki, and K. Aihara: Identifying pre-disease signals before metabolic syndrome in mice by dynamical network biomarkers, Scientific Reports, 9:8767 (2019). https://doi.org/10.1038/s41598-019-45119-w

```python
# load data
def load_data():
  file_name = '../data/data.tsv'
  data_df = pd.read_csv(file_name, sep='\t', header=[0,1,2], index_col=0)
  data_df.columns = data_df.columns.droplevel(0)
  data_df = (data_df / stats.trim_mean(data_df, 0.02)).apply(np.log2)
  return data_df

def load_info():
  file_name = '../data/info.tsv'
  return pd.read_csv(file_name, sep='\t')

if __name__ == '__main__':
  data_df = load_data()
  info_df = load_info()
    
# process data
import numpy as np
import pandas as pd

def preprocess_data():

  # NOTE: the mapping relations between probe names and gene symbols
  # were extracted from 028005_D_GeneList_20171030.txt obtained at:
  # https://earray.chem.agilent.com/earray/catalogGeneLists.do?action=displaylist
  # It was the latest at the momemnt of the paper submission.
  data_file_name = '../data/Your_Matrix_data.txt'
  mapping_file_name = '../data/id_meta_mapping.tsv' 
  output_file_name = '../data/data.tsv'

  # load expression data file ---------------------------------

  df = pd.read_csv(data_file_name, sep='\t', index_col=0, 
                   header=[35,71]) # '!Sample_title', 'ID_REF'
  df.drop('!series_matrix_table_end', inplace=True) # => 62976 x 64
  df.index = df.index.astype(int)

  # rename columns ------------------------------------------------

  col_idx         = df.columns.get_level_values(0)
  condition_idx   = col_idx.str.split('_').str[0]
  week_idx        = col_idx.str.split('_').str[1].str.replace('w','')
  df.columns = [np.arange(len(col_idx)), condition_idx, week_idx]
  df.columns.names = ['sample_no', 'condition', 'week']
  
  # id conversion -------------------------------------------------

  mapping_sr = pd.read_csv(mapping_file_name, sep='\t',
                           usecols=['FeatureNum', 'GeneSymbol'],
                           index_col=0, squeeze=True).dropna()
  df = df.loc[mapping_sr.index].rename(mapping_sr)
  df = df.groupby(axis=0,level=0).mean()
  df.index.name = 'gene_symbol'

  # save to file -------------------------------------------------

  df.to_csv(output_file_name, sep='\t')


if __name__ == '__main__':
  preprocess_data()
```

We use the jupyter notebook to run the script, note that you need to put the Matrix-data and Meta data in the same directory.

## Ensembl id & Symbol id

Sometimes, the expression matrix always has the row names, however those are Ensembl id......

We can't say that Ensembl id is bad, but in some way Symbol id is indeed a better choice. So I have to transform Ensembl id into Symbol id.

I have to clarify that both Ensembl & Symbol id are all Gene id. They are completely different from Probe id used in Chip data. Let's make it simble:

We can use Ensembl id or Symbol in expression matrix. But we can never use probe id in a normal expression matrix.

In order to make this transformation, I have to seek help in Ensembl main website, and I find something called Biomart.

The steps are listed below:

```shell
# 1.Extract Ensembl ID from Expression(Not in Excel!!!) 
cat Your_Expr_Matrix | awk '{print $1}' | sed -i '1d'> EnsemblID.txt

-----------------------------------------------------------------------------------------------------

# Or use R
rt<-read.table(file = 'Your_Expr_Matrix.txt',header = T,row.names = 1)
EnsemblID<-row.names(rt)
write.csv(EnsemblID,file = 'EnsemblID.txt')
# Next though we can use biomart API in R, but it's easier to use website directely
```

We use biomart to retrieve Gene Symbol matched the Ensembl id:

[Biomart](http://grch37.ensembl.org/biomart/martview/ebd624a5c1df8c3f5c571a8afb446513)

And  then I upload my gene list, and I can easily get the gene feature whatever I wanted. Together with Function merge in R, we can then transform Ensembl id into Symbol.

However, this is not a general way to solve this problem, I will then get a general method using R and biomart API.

