ReTimeML is a bespoke, freeware, application for estimating the retention times (RTs) of
ceramide and sphingomyelin sphingolipids from complex total ion chromatograms. ReTimeML
works on the foundation that every LC-MS/MS experiment will have pre-determined or known
RTs as a consequence of the internal controls, external calibrators or quality control mixtures
used in their analysis. Employing these RTs as reference points, ReTimeML can extrapolate
unknowns using its machine-learned regression library of mass-to-charge (m/z) vs relative RTs
for ceramides and sphingomyelins.
ReTimeML provides users with a 2-D profile of its estimations (as a plot of m/z vs estimated RT),
together with a list for all sphingolipids RTs listed, which can be downloaded as a .csv file or
directly copied into Excel or a similar spreadsheet.

Basics
To operate you will need to upload a separate data file for ceramides and sphingomyelins of
interest. Input should have &quot;Train&quot; and &quot;Test&quot; labelled in a .csv as below.
-file is in csv format -file has a column &quot;Lipid ID&quot; (without the quotations) which specifies all the lipid
names in LipidMaps TM nomenclature (e.g. Cer (d18:1/20:1)). -file has a column &quot;Retention Time&quot;
(without the quotations) which specifies the user-recorded retention times -file has a column &quot;Type&quot;
(without the quotations) with values either &quot;Train&quot; or &quot;Test&quot; (without the quotations) - note any -record
with &quot;Train&quot; must have a corresponding retention time as this refers to a known value that will be used
to determine the unknowns (“Test”)
Uploaded data (drag and drop option) triggers the automatic calculation of RTs. Run your
sphingolipid list and view the predicted retention times you are interested in.

*We are currently in the process of expanding ReTimeMLs capacity to other lipids and welcome support from
colleagues and experts in the field willing to share their prior datasets of verified lipids, with corresponding
precursor mass and RT, which can be uploaded via the Google form: https://forms.gle/p9fGuzuujZqDhcpD6.
University of Sydney (Created by Michael Allwright and Timothy Couttas).