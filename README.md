# wingflow_grasshoppers
Matlab code for analyzing particle trajectories from fluoresecent data set. 

Code doesn't have a particular structure, but does sort particles into a normalized wing region and into the five wing regions listed in the paper: leading edge, membrane, lattice, wing tip, trailing edge.

Within the data_sheets.xlsx are several tables: 

R - a matrix which initially sorts the data:
  Col1: k
  Col2: maximum velocities
  Col3 and Col4: nothing
  Col5: alpha_1 (the first part of the Womersley equation, 2nd part calculated in the wing-region part of the code)
  Col6: left (noted by 0) or right wings (noted by 1)
  Col7: particle mean
  Col8: particle max
  Col9: particle median velocity
  Col10: re_1 (first part of Reynolds number equation - 2nd part calculated in the wing-region part of the code)
  Col11: cum_dist (cumulative distance a particle travels)
  Col12: numDig (number of frames digitized per trajectory)
  Col13: numPks (number of peaks within a velocity trace)
  Col14: time int (time interval per sec per image slice)
  Col15: regiondata (which region a given particle is in)
  
 The other tables in the data sheet are: 
 wingReg - particle velocity means
 wingRegMed - particle velocity medians
 wingRegMax - particle velocity maximums
 wingRegFreq - particle frequencies per region
 wingRegWo -  Womersley number per region
 wingRegRE -  Reynolds number per region
 wingRegPe - Peclet number per region
 
 All of these tables have 10 columns and are, FW Leading Edge,	FW Membrane,	FW Wing Tip,	FW Lattice, FW Trailing Edge,	HW Leading Edge,	HW Membrane,	HW Wing Tip,	HW Lattice, and HW Trailing Edge (where FW - forewing and HW - hindwing). 


The last two tables - regions and regionsMedian - have 28 columns
Category order:'ABD','PRO','DH','SCUT','FW Hinge','FW','HW Hinge','HW'

These columns were used to calculate velocity means and medians, and to parse out specific locations in these body regions:
Col1: FW Hinge mean velocity
Col2: X-location of FW Hinge
Col3: Y-location of FW Hinge
Col4: HW Hinge mean velocity
Col5: X-location of HW Hinge
Col6: Y-location of HW Hinge
Col7: Abdomen mean velocity
Col8: X-location of Abdomen
Col9: Y-location of Abdomen
Col10: Dorsal Heart mean velocity
Col11: X-location of Dorsal Heart
Col12: Y-location of Dorsal Heart
Col13: Scutellar Branch mean velocity
Col14: X-location of Scutellar Branch
Col15: Y-location of Scutellar Branch
Col16: Pronotum mean velocity
Col17: X-location of Pronotum
Col18: Y-location of Pronotum
Col19: Forewing mean velocity
Col20: nothing (variable was unused)
Col21: nothing (variable was unused)
Col22: X-location of Forewing
Col23: Y-location of Forewing
Col24: Hindwing mean velocity
Col25: nothing (variable was unused)
Col26: nothing (variable was unused)
Col27: X-location of Hindwing
Col28: Y-location of Hindwing



