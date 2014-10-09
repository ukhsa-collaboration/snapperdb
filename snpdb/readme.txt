This module will take

* config file with reference genome, snpdb name, pg host, pg uname, pg pword
* path to vcf (assuming that within our pipeline, depth has been checked before vcf made, within other pipelines, that depth
has been check elsewhere - will also go into strain stats (calculate from vcf) so can select based on this in query anyway).


This module will produce

* log of output
* new db entry