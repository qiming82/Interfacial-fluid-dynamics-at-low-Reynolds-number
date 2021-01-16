#FC=pgf90
FC=gfortran
#FFLAGS=-g
# capillary breakup
#---------
#
OBJ = arc_evl.o bess_I01K01.o df.o dgesv.o ell_int.o ell_int2.o gauss_leg.o gauss_log.o \
gsint.o drop_sfct_HAS.o polyn.o regridf.o regridf_adp.o regridf_sft2.o regridf_sft_adp.o getsfct0.o getsfct2.o\
spline.o seval.o sgf_ax_fs.o sgf_ax_g.o sgf_ax_str.o splc_clm.o strint.o strint2.o thomas.o
#
FFLAGS = -O3 -g
#
break: $(OBJ)
	$(FC) $(FFLAGS) -o break $(OBJ)	
drop_sfct_HAS.o: drop_sfct_HAS.for
	$(FC) $(FFLAGS) -c drop_sfct_HAS.for
arc_evl.o: arc_evl.for
	$(FC) $(FFLAGS) -c arc_evl.for
bess_I01K01.o: bess_I01K01.f
	$(FC) $(FFLAGS) -c bess_I01K01.f
df.o: df.for
	$(FC) $(FFLAGS) -c df.for
#initl.o: initl.for
#	$(FC) $(FFLAGS) -c initl.for
dgesv.o: dgesv.f
	$(FC) $(FFLAGS) -c dgesv.f
ell_int.o: ell_int.f
	$(FC) $(FFLAGS) -c ell_int.f
ell_int2.o: ell_int2.for
	$(FC) $(FFLAGS) -c ell_int2.for
gauss_leg.o: gauss_leg.f
	$(FC) $(FFLAGS) -c gauss_leg.f
gauss_log.o: gauss_log.for
	$(FC) $(FFLAGS) -c gauss_log.for
gsint.o: gsint.for
	$(FC) $(FFLAGS) -c gsint.for
getsfct0.o: getsfct0.for
	$(FC) $(FFLAGS) -c getsfct0.for
#getsfct0.o: getsfct0.for
#	$(FC) $(FFLAGS) -c getsfct0.for
getsfct2.o: getsfct2.for
	$(FC) $(FFLAGS) -c getsfct2.for
polyn.o: polyn.for
	$(FC) $(FFLAGS) -c polyn.for
regridf.o: regridf.for
	$(FC) $(FFLAGS) -c regridf.for
regridf_adp.o: regridf_adp.for
	$(FC) $(FFLAGS) -c regridf_adp.for
#regridf_ax0.o: regridf_ax0.for
#	$(FC) $(FFLAGS) -c regridf_ax0.for
#regridf_ax1.o: regridf_ax1.for
#	$(FC) $(FFLAGS) -c regridf_ax1.for
regridf_sft2.o: regridf_sft2.for
	$(FC) $(FFLAGS) -c regridf_sft2.for
regridf_sft_adp.o: regridf_sft_adp.for
	$(FC) $(FFLAGS) -c regridf_sft_adp.for
#regridf_sft_adpb.o: regridf_sft_adpb.for
#	$(FC) $(FFLAGS) -c regridf_sft_adpb.for
spline.o: spline.f
	$(FC) $(FFLAGS) -c spline.f
seval.o: seval.f
	$(FC) $(FFLAGS) -c seval.f
sgf_ax_fs.o: sgf_ax_fs.f
	$(FC) $(FFLAGS) -c sgf_ax_fs.f
sgf_ax_g.o: sgf_ax_g.for
	$(FC) $(FFLAGS) -c sgf_ax_g.for
sgf_ax_str.o: sgf_ax_str.for
	$(FC) $(FFLAGS) -c sgf_ax_str.for
splc_clm.o: splc_clm.f
	$(FC) $(FFLAGS) -c splc_clm.f
strint.o: strint.for
	$(FC) $(FFLAGS) -c strint.for
strint2.o: strint2.for
	$(FC) $(FFLAGS) -c strint2.for
thomas.o: thomas.f
	$(FC) $(FFLAGS) -c thomas.f
#
#clean:
#\rm -rf *.o pinch
#cleaner: clean
#\rm -rf output/*
#
