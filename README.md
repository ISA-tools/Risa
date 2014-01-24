<p align="center">
<img src="http://isatab.sourceforge.net/assets/img/tools/tools-table-images/risatab.png" align="center" alt="Risa"/>
</p>


Risa is an R package that is part of the ISA tools suite (http://isa-tools.org). Risa supports parsing, saving and updating ISA-tab datasets. It also builds bridges from the ISA-Tab syntax to analysis pipelines for specific assay types, such as mass spectrometry and DNA microarray assays, by building R objects from the metadata required for other packages downstream, such as xcms and affy, respectively. In addition, Risa includes functionality to suggest packages in BioConductor that might be relevant for the assay types in the ISA-TAB dataset being considered. This recommentation functionality relies on the BioCViews annotations provided by each BioConductor package.



The Risa package is available in Bioconductor:
  - [release version](http://www.bioconductor.org/packages/release/bioc/html/Risa.html)
  - [development version](http://www.bioconductor.org/packages/devel/bioc/html/Risa.html)

For more information about the ISA tools, consider: 

- General info: <http://isa-tools.org>
- Tools' overview in this short paper: <http://bioinformatics.oxfordjournals.org/content/26/18/2354.full.pdf+html>
- Issue tracking and bug reporting: <https://github.com/ISA-tools/Risa/issues>
- Mainline source code: <https://github.com/ISA-tools/Risa>
- Twitter: [@isatools](http://twitter.com/isatools)
- IRC: [irc://irc.freenode.net/#isatab](irc://irc.freenode.net/#isatab)
- [Development blog](http://isatools.wordpress.com) 

### Read the Paper in BMC Bioinformatics!
Access the Open Access BMC Bioinformatics article on Risa [here](http://www.biomedcentral.com/1471-2105/15/S1/S11).

Alejandra González-Beltrán, Steffen Neumann, Eamonn Maguire, Susanna-Assunta Sansone and Philippe Rocca-Serra.  
The Risa R/Bioconductor package: integrative data analysis from experimental metadata and back again 
BMC Bioinformatics 2014, 15(Suppl 1):S11  [doi:10.1186/1471-2105-15-S1-S11](http://dx.doi.org/10.1186/1471-2105-15-S1-S11)


## Development

If you have feature requests or find any issues when using the tracker, please let us know through the issues tracker at [https://github.com/ISA-tools/Risa/issues]. 

### Contributing

You should read this article about Git Flow: <http://scottchacon.com/2011/08/31/github-flow.html>. It's a really useful tutorial on how to use Git for collaborative development.

1. Fork it.
2. Clone your forked repository to your machine
3. Create a branch (`git checkout -b myRisa`)
4. Make your changes
5. Run the tests (`mvn clean test`)
6. Commit your changes (`git commit -am "Added something useful"`)
7. Push to the branch (`git push origin myRisa`)
8. Create a [Pull Request](http://help.github.com/pull-requests/) from your branch.
9. Promote it. Get others to drop in and +1 it.

