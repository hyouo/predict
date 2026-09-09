# v7 comparator strengthening before new scoring

No new outer test scores have been computed. Alongside inherited width=1 RBF, include an RBF bandwidth grid {0.25,1,4}. The same bandwidth grid applies to every target/pathway and null kernel mixture, so a prior cannot receive a bandwidth advantage. The full linear comparator remains. Ridge grid and analytic [0,1] amplitude are unchanged. Report the larger eta search space of prior families.

Exact-ID mapping audit has 904 unambiguous target-annotated drugs among 1,726, 128 matched names without targets, 688 unmatched CIDs, and 6 incompatible-CID annotations excluded. There are 971 distinct accepted target symbols, of which 108 are measured landmarks. Direct priors cover 276 drugs; the prespecified human Reactome expansion covers 814. Source/sample/pathway hashes are in the first protocol. Degree-preserving random graph seed 1707 completed 62,328 accepted swaps among 3,337 edges, preserving each drug's measured/unmeasured target counts and each target's frequency.

The residual target-operator uses rank<=6 source baseline PCs and rank<=16 normalized target incidence singular vectors. RBF residual labels come from leaving each fitted source context out, including re-estimating baseline scale. Main RBF kernel choices are selected only on the fixed validation contexts. Target-operator correction strength zero is included. A shuffled-target operator has the same ranks and searches. Predictors with no direct/pathway information explicitly fall back to the matching ordinary kernel; unknown drug annotation is never treated as a zero drug response.

This is explicitly retrospective development on exposed L1000 contexts. Neither these commits nor prediction hashes restore first-use blind status.
