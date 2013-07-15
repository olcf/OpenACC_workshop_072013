mg_v00_orig.c

This is the code as downloaded with ragged 3d and 4d arrays.

mg_v06_orig.c

This is a full acceleration using OpenACC, where all the arrays are manually
deep-created on the GPU. This is really slow, so only run this code compiled
with CLASS=S or CLASS=W.


mg_v00.c

This rewrites mg_v00_orig.c so that the 3d data arrays are contiguous blocks
of memory, with the array indexing flattened to one index.

The 4d arrays (u,r) are ragged on the outer, level index but the 3d data arrays that
are indexed by this are made continuous.

mg_v0?.c

Develops the OpenACC port incrementally

mg_v0?a.c

Are tuning variations, sorting out minor scheduling issues.
