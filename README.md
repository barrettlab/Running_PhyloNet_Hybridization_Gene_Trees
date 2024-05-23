# Part I. Running_PhyloNet_Hybridization_Gene_Trees

Preparing gene tree file and running Phylonet in R

See the tutorial and citations therein (https://phylogenomics.rice.edu/html/phylonetTutorial.html)

 0. Load libraries in R
```{r}
library(phytools)
```
 1. in R, read in gene trees as a multiphylo object
```{r}
tre <- read.tree("treeshrink_bestparalog.treefile")
```
 2. check to see that it is a multiphylo object
```{r}
class(tre)

"multiphylo"
```

 3. specify tips to be pruned from tree
```{r}
taxa_to_drop<_c("fastp_Coelia_triptera_S10","fastp_Yoania_prainii_S14","fastp_Yoania_japonica_S12","fastp_Calypso_bulbosa_var_americana_S4","fastp_Calypso_bulbosa_Asia_S37","fastp_Changnienia_amoena_S44","fastp_Tipularia_japonica_S45","fastp_Corallorhiza_maculata_var_maculata2_S36","fastp_Corallorhiza_mertensiana_S30","fastp_Corallorhiza_maculata_var_mexicana_S17","fastp_Corallorhiza_macrantha_S21","fastp_Corallorhiza_bulbosa_S19","fastp_Corallorhiza_odontorhiza_var_odontorhiza2_S38","fastp_Corallorhiza_odontorhiza_var_pringlei_S25","fastp_Corallorhiza_odontorhiza_Mexico_S47","fastp_Corallorhiza_wisteriana_Eastern_US_S27","fastp_Corallorhiza_striata_var_striata_S5","fastp_Corallorhiza_striata_Sierra_Nevada_S9","fastp_Corallorhiza_striata_CA_Coast_Ranges_S11","fastp_Corallorhiza_involuta_S1","fastp_Oreorchis_indica2_S46","fastp_Cremastra_variabilis_S24","fastp_Cremastra_aphylla_S18","fastp_Cremastra_saprophytica_S22","fastp_Cremastra_unguiculata_S43","fastp_Aplectrum_hyemale_S8","fastp_Govenia_superba_S6","fastp_Govenia_capitata_S42","fastp_Dactylostalix_ringens_S39","fastp_Ephippianthus_schmidtii_S26","fastp_Ephippianthus_sawadanus_S31","fastp_Brassavola_glauca_S40")
```
 4. drop all tips except the ones you want to keep
```{r}
pruned.tree<-drop.tip.multiPhylo(tre,taxa_to_drop)
```
 5. plot the first gene tree to see if it worked
```{r}
plot(pruned.tree[[1]])
```
 6. write the tree to file
```{r}
write.tree(pruned.tree,file="treeshrink_Calyps_subset.nex")
```
 7. You need to convert the tree to nexus format, either manually or via phytools:
```{r}
writeNexus(pruned.tree,file="nexustest.nex")
```
 8. You'll need to edit the nexus file to look like this:

```bash

#NEXUS
 
BEGIN TREES;


Tree gt001=(((fastp_Oreorchis_coreana_S32:0.0176422487,fastp_Cremastra_appendiculata_S20:0.0085695487):1e-06,(((fastp_Oreorchis_indica1_S35:2e-06,fastp_Oreorchis_fargesii_S34:1e-06):1e-06,fastp_Oreorchis_bilamellata_S48:1e-06):2.0377e-06,(fastp_Oreorchis_erythrochrysea_S33:2e-06,fastp_Oreorchis_patens_S29:1e-06):1e-06):2.0306e-06):1e-06,(((fastp_Corallorhiza_bentleyi_S3:2.5214e-06,(fastp_Corallorhiza_striata_var_vreelandii_S7:3.0243e-06,fastp_Corallorhiza_trifida_S2:1e-06):1e-06):0.008596816,fastp_Corallorhiza_wisteriana_Western_US_S28:0.0085942908):2.0473e-06,fastp_Corallorhiza_maculata_var_occidentalis_S15:2e-06):0.0086608932);
Tree gt002=(((((fastp_Corallorhiza_bentleyi_S3:0.0037893074,((fastp_Corallorhiza_wisteriana_Western_US_S28:0.0096400708,fastp_Corallorhiza_maculata_var_occidentalis_S15:0.0018930636):0.0018905778,(fastp_Corallorhiza_striata_var_vreelandii_S7:0.0018965563,fastp_Corallorhiza_trifida_S2:0.0076179904):0.0018905325):1.001e-06):0.001890477,fastp_Oreorchis_coreana_S32:0.0038012431):1.001e-06,((fastp_Oreorchis_erythrochrysea_S33:0.0038074333,fastp_Oreorchis_indica1_S35:2.002e-06):0.0018910937,fastp_Oreorchis_patens_S29:0.0037836267):1.001e-06):1.001e-06,(fastp_Oreorchis_bilamellata_S48:1.001e-06,fastp_Oreorchis_fargesii_S34:0.0018865137):0.0018889098):0.0037844083,fastp_Cremastra_appendiculata_S20:0.0037834846);
Tree gt003=((((fastp_Oreorchis_bilamellata_S48:1.1108e-06,fastp_Oreorchis_coreana_S32:1.1108e-06):0.0041928447,fastp_Oreorchis_fargesii_S34:1.1108e-06):0.0190836901,(fastp_Corallorhiza_bentleyi_S3:0.0305309322,(fastp_Corallorhiza_trifida_S2:0.0126344377,(fastp_Oreorchis_erythrochrysea_S33:0.0045712793,fastp_Oreorchis_indica1_S35:0.0083840509):0.0041861576):1.1108e-06):0.0049402476):0.0167464358,fastp_Cremastra_appendiculata_S20:0.0219480731);
Tree gt004=(fastp_Cremastra_appendiculata_S20:0.0030859065,((((fastp_Corallorhiza_bentleyi_S3:1.0004e-06,fastp_Corallorhiza_striata_var_vreelandii_S7:0.0060026864):0.0029992081,fastp_Corallorhiza_trifida_S2:0.0069111331):0.0019470118,(((fastp_Oreorchis_bilamellata_S48:1.0004e-06,fastp_Oreorchis_fargesii_S34:0.0034964139):0.0035081263,fastp_Oreorchis_coreana_S32:0.0117322654):1.0004e-06,fastp_Oreorchis_patens_S29:1.0004e-06):0.0038676149):1.0004e-06,fastp_Oreorchis_indica1_S35:0.003441931):0.002733112);
Tree gt005=(((((fastp_Cremastra_appendiculata_S20:0.0146515673,((fastp_Corallorhiza_bentleyi_S3:0.00368338,(fastp_Oreorchis_bilamellata_S48:0.0057207567,((fastp_Oreorchis_coreana_S32:1.0322e-06,fastp_Oreorchis_fargesii_S34:1.0322e-06):1.0322e-06,fastp_Oreorchis_patens_S29:1.0322e-06):1.0322e-06):0.0074283596):1.0322e-06,(fastp_Corallorhiza_trifida_S2:1.0322e-06,fastp_Oreorchis_erythrochrysea_S33:1.0322e-06):2.0644e-06):1.0322e-06):0.0035918253,fastp_Oreorchis_indica1_S35:1.0322e-06):3.4491e-06,fastp_Corallorhiza_maculata_var_occidentalis_S15:1.0322e-06):1.0322e-06,fastp_Corallorhiza_wisteriana_Western_US_S28:1.0322e-06):3.7869e-06,fastp_Corallorhiza_striata_var_vreelandii_S7:1.0322e-06);

...the rest of the gene trees...

END;
 
BEGIN PHYLONET;
 
InferNetwork_MPL (all) 0 -pl 30;
 
END;

```

 9. In the above, you need the PHYLONET block, with this general format

```bash

END;
 
BEGIN PHYLONET;

Algorithm (taxa_to_analyze) #hybridizations -pl <#threads>

END;

```

 10. Now, you are ready to run PHYLONET. These analyses for each H-value (# hybridizations) take 10-30 minutes.
 11. Running the block above specifies zero hybridizations, and essentially finds the "species tree"
```bash
java -jar /usr/local/src/PhyloNet.jar treeshrink_Calyps_subset.nex
```
 12. When finished, you'll get a bunch of output. Take the first of the final five trees and the likelihood score to calculate the AIC
```bash

Inferred Network #1:
((((((fastp_Oreorchis_fargesii_S34:1.0,fastp_Oreorchis_bilamellata_S48:1.0):0.27050983124842287,(fastp_Oreorchis_patens_S29:1.0)#H1:1.0::0.981174812035372):0.026811468074220527,fastp_Oreorchis_coreana_S32:1.0):0.3289465597576068,(((((#H1:1.0::0.018825187964627954,fastp_Corallorhiza_striata_var_vreelandii_S7:1.0):5.910541841982505,fastp_Corallorhiza_bentleyi_S3:1.0):0.33548821973179177,(fastp_Corallorhiza_trifida_S2:1.0,(fastp_Corallorhiza_wisteriana_Western_US_S28:1.0,fastp_Corallorhiza_maculata_var_occidentalis_S15:1.0):0.32911936196475744):0.1732033524782276):0.47485874357798835,(fastp_Oreorchis_indica1_S35:1.0,fastp_Oreorchis_erythrochrysea_S33:1.0):0.16451242730249852):0.2193859870956898)#H2:0.14515715991898692::0.7976364616620312):0.7588902614642841,fastp_Cremastra_appendiculata_S20:1.0):0.7544223398620533,#H2:0.7338327544302619::0.20236353833796883);
Total log probability: -35154.69555932267
Visualize in Dendroscope : ((((((fastp_Oreorchis_fargesii_S34,fastp_Oreorchis_bilamellata_S48),(fastp_Oreorchis_patens_S29)#H1),fastp_Oreorchis_coreana_S32),(((((#H1,fastp_Corallorhiza_striata_var_vreelandii_S7),fastp_Corallorhiza_bentleyi_S3),(fastp_Corallorhiza_trifida_S2,(fastp_Corallorhiza_wisteriana_Western_US_S28,fastp_Corallorhiza_maculata_var_occidentalis_S15))),(fastp_Oreorchis_indica1_S35,fastp_Oreorchis_erythrochrysea_S33)))#H2),fastp_Cremastra_appendiculata_S20),#H2);
```

 13. Repeat the analysis for however many numbers of H (0,1,2,3,4,5,...). You need to directly edit the nexus file to do this.

 14. Save the likelihood scores as a column in excel (or do this in R).

 15. The number of parameters (k) = the number of total internal + terminal branches in the tree PLUS the number of pre-specified hybridization events (k + H). Use these for AIC calcs, where AIC = 2*k - 2*lnL.

 16. Calulate AIC scores in excel with " =((2*C3)-(2*LN(B3)))" where the likelihood is in cell B3 and k is in B4.

 17. Calculate the delta AIC (AIC for value of H - minAIC) in excel.

 18. Calculate AIC weights (wAIC) as "=EXP(-0.5*E3)" where the deltaAIC is in cell E3, and drag down.

### Voila! 

### The H-values (# of hybridizations) with the lowest AIC and highest wAIC is the optimal. It could very well be zero.

# Part II. Analysis with PhyloNetworks

Following (https://crsl4.github.io/PhyloNetworks.jl/dev/man/snaq_plot/)

 1. Start julia and add phylonetworks
```bash
julia
using PhyloNetworks
```
 2. Read in gene trees, view tree #3
```bash
genetrees = readMultiTopology("calyps_phylonetworks.tre");
genetrees[3]
```
 3. Load phyloplots and plot tree #3
```bash
using PhyloPlots
plot(genetrees[3]); # tree for 3rd gene
```
 4. Calculate quartet Concordance Factors
```bash
q,t = countquartetsintrees(genetrees);
```
 5. Write CFs to a table to save and view
```bash
using CSV
df = writeTableCF(q,t)   # data frame with observed CFs: gene frequencies
CSV.write("tableCF.csv", df); # to save the data frame to a file
raxmlCF = readTableCF("tableCF.csv") # read in the file and produces a "DataCF" object
less("tableCF.csv")

raxmlCF = readTrees2CF(genetrees, whichQ="rand", numQ=200, CFfile="tableCF10.txt")
```
 6. Get a starting tree -- in this case, the astral "species tree"
```bash
astraltree = readTopology("astral.tre")
```
 7. For multithreading/parallel jobs
```bash
using Distributed
addprocs(30)
@everywhere using PhyloNetworks
```
 8. Run the first network analysis with H=0, essentially find the "species tree"
```bash
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234)
```
 9. Got booted out of julia for some reason, needed to start over!
```bash
using PhyloPlots
q,t = countquartetsintrees(genetrees); # read in trees, calculate quartet CFs
df = writeTableCF(q,t)

using CSV
CSV.write("tableCF.csv", df);
raxmlCF = readTableCF("tableCF.csv")
raxmlCF = readTrees2CF(genetrees, whichQ="rand", numQ=200, CFfile="tableCF10.txt")

astraltree = readTopology("astral.tre")
plot(astraltree, showedgelength=true);
```
 10 Run the first network analysis with H=0, essentially find the "species tree"
```bash
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234)
```
 11. Copy the log-likelihood from the screen output after run finishes, or get from "net0.out" file

 12. Run the second network analysis with H=1, allowing 1 hybridization event, using 'net0' as starting tree
```bash
net1 = snaq!(net0, raxmlCF, hmax=1, filename="net1", seed=2345) # this runs for ~30 min
```
 13. Check out the output
```bash
plot(net1, showgamma=true);
less("net1.err") # would provide info about errors, if any
less("net1.out") # main output file with the estimated network from each run
less("net1.networks") # extra info
net1
```
 14. Save the network file
```bash
writeTopology(net1, round=true, digits=2)
```
 15. Run the third network analysis with H=2, allowing 2 hybridization events, using 'net0' as starting tree
```bash
net2 = snaq!(net0,raxmlCF, hmax=2, filename="net2", seed=3456)
plot(net2, showgamma=true);
```
 16. Run the rest of the analyses analysis with H=3-5, allowing 3-5 hybridization events, using 'net0' as starting tree
### Save the likelihoods for AIC calcs

```bash
net3 = snaq!(net0,raxmlCF, hmax=3, filename="net3", seed=4567)
net4 = snaq!(net0,raxmlCF, hmax=4, filename="net4", seed=1437)
net5 = snaq!(net0,raxmlCF, hmax=5, filename="net5", seed=8701)
```




