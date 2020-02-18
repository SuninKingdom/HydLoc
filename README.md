<h1 align="center">
  <img alt="logo" src="logo.png">
</h1>
<h1 align="center">
  <img alt="internal title" src="internal title.svg">
</h1>


## Table of Contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Notice](#notice)
* [Citation](#citation)

### Introduction:

HydLoc is a random forest based hydroxylation sites predictor for human proteins. HydLoc is provided as a 
Microsoft Windows executable packed up with PyInstaller, which makes it easy to use. The source code of HydLoc is also uploaded. 
It is a python script called 'HydLoc.py'. Users can run it to predict hydroxylation sites too. Because the training dataset of HydLoc 
is human proteins, it is recommended that HydLoc is only used to predict human hydroxylation sites.

### Installation:

Download the package. The size of 'HydLoc.exe' is about **283 MB**. 
It can not be downloaded directly through 'HydLoc-master.zip'. 
Please download the complete file in the “HydLoc.exe” item on GitHub 
(https://github.com/SuninKingdom/HydLoc/blob/master/HydLoc.exe). 
Make sure that **'clf1', 'clf2', 'HL.ico'** are in the same folder with **'HydLoc.exe'**, 
and do not change the name of them. For the beauty and stability of the interface, the version 
of Office is better not earlier than 2016.

### Usage:

* 1.Start HydLoc.exe by double click.

* 2.Enter or paste the query protein sequences in the text box, or upload a file with the query sequences.

* 3.Click the 'Submit' button.

### Notice:

* 1.The protein sequence should be given in Fasta format.

* 2.If there exists amino acids beyond the 20 standard amino acids in the query sequences, HydLoc will replace 
    them with native amino acids picked up randomly. 

* 3.If the proline/lysine residues appear in the first or the last nine/fourteen in the query sequence, 
    they will be regarded as non-hydroxylated sites by HydLoc. 

* 4.If the protein query list was too long, more patience is needed for the result.And if the result could 
    not be completely projected on the result canvas, please 'Save as' the result, the whole outcome may in
    the saved file.

### Citation:

Please cite the following paper for using HydLoc:
Qixing Huang, Xingyu Chen, Yang Wang, Jinlong Li, Haiyan Liu,Yun Xie, Zong Dai, Xiaoyong Zou* and Zhanchao Li*. 
HydLoc: a tool for hydroxyproline and hydroxylysine sites prediction in the human proteome. *CHEMOMETR INTELL LAB*, 2019.
