nohup root -b -l -q rootlogon.C 'FillTempHist.cpp("Be")' > & TempBewide_L1Inner.log &
nohup root -b -l -q rootlogon.C 'FillTempHist.cpp("B")' > & bkgTempB_frag_L1Innertest.log &

rootfile log
_temp_wide_bkg_final, 25.5.15, iss, L1L2 valid study, bkg selection, richq>3 for mass dis., z/a = 5/10 cutoff beta
_temp_wide_bkg_frag,  25.5.15, mc,  L1L2 valid study, bkg selection, richq>3 for mass dis., L1 Boron(sel. by truth[0] or cut[1]) frag to Inner Be using mtrpar[0],[1]
_temp_wide_bkgcut,    25.5.15, mc,  be mc mass temp with bkgcut(no tof and l1 q, inner 3.5~4.5, richq>3, unbiasedL1, and other)

_temp_wide_bkg_frag_fitstudy 25.5.20 same as _temp_wide_bkg_frag, add hist to record id in each ek, and pure frag be7 9 10 mass dis.
*_4.8 25.5.20 unbiasedL1 q 4.8-5.5
*_4.8to5.4 25.5.20 unbiasedL1 q 4.8-5.4
*_SDIATQ 25.5.20 unbiasedL1 q 5.0-5.4

above all inner-unbiasedL1 acc

"_temp_wide_L1Inner.root" : "_temp_wide_bkg_frag_L1Inner.root" 25.6.9, For Be, Mass Hist with All L1Inner + Beta Det Cuts. For Bor, Bkg 4.8-5.5 L1Inner Sample Study