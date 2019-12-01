
for(int i=0;i<NGHOSTS;i++) { free(bcstruct.outer[i]);  free(bcstruct.inner[i]); }
free(bcstruct.num_ob_gz_pts); free(bcstruct.num_ib_gz_pts);
