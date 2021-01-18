FC=pgf90
#FFLAGS=-g
# caff
#---------
#
OBJ = caf_main.o arc_evl.o bess_I01K01.o df.o dgesv.o ell_int.o ell_int2.o gauss_leg.o gauss_log.o \
gsint_ct.o gsint.o lgf_ax_fs.o lgf_ax_gp_ct.o lgf_ax_gpc_ct.o polyn.o regridf.o regridf2.o \
seval.o sgf_ax_fs.o sgf_ax_1p.o sgf_ax_1p_ct.o splc_clm.o thomas.o \
int_pot_s.o int_pot_sc.o geten.o strint.o
#
FFLAGS = -O3 -g
#
CAfflow: $(OBJ)
	$(FC) $(FFLAGS) -o CAfflow $(OBJ)	
caf_main.o: caf_main.for
	$(FC) $(FFLAGS) -c caf_main.for
arc_evl.o: arc_evl.for
	$(FC) $(FFLAGS) -c arc_evl.for
bess_I01K01.o: bess_I01K01.f
	$(FC) $(FFLAGS) -c bess_I01K01.f
df.o: df.for
	$(FC) $(FFLAGS) -c df.for
dgesv.o: dgesv.f
	$(FC) $(FFLAGS) -c dgesv.f
#dgmres.o: dgmres.f
#	$(FC) $(FFLAGS) -c dgmres.f
ell_int.o: ell_int.f
	$(FC) $(FFLAGS) -c ell_int.f
ell_int2.o: ell_int2.for
	$(FC) $(FFLAGS) -c ell_int2.for
gauss_leg.o: gauss_leg.f
	$(FC) $(FFLAGS) -c gauss_leg.f
gauss_log.o: gauss_log.for
	$(FC) $(FFLAGS) -c gauss_log.for
gsint_ct.o: gsint_ct.for
	$(FC) $(FFLAGS) -c gsint_ct.for
gsint.o: gsint.for
	$(FC) $(FFLAGS) -c gsint.for
#guass_el.o: guass_el.for
#	$(FC) $(FFLAGS) -c guass_el.for
lgf_ax_fs.o: lgf_ax_fs.f
	$(FC) $(FFLAGS) -c lgf_ax_fs.f
lgf_ax_gp_ct.o: lgf_ax_gp_ct.for
	$(FC) $(FFLAGS) -c lgf_ax_gp_ct.for
lgf_ax_gpc_ct.o: lgf_ax_gpc_ct.for
	$(FC) $(FFLAGS) -c lgf_ax_gpc_ct.for
int_pot_s.o: int_pot_s.for
	$(FC) $(FFLAGS) -c int_pot_s.for
int_pot_sc.o: int_pot_sc.for
	$(FC) $(FFLAGS) -c int_pot_sc.for
geten.o: geten.for
	$(FC) $(FFLAGS) -c geten.for
#mat_it.o: mat_it.for
#	$(FC) $(FFLAGS) -c mat_it.for
polyn.o: polyn.for
	$(FC) $(FFLAGS) -c polyn.for
regridf.o: regridf.for
	$(FC) $(FFLAGS) -c regridf.for
regridf2.o: regridf2.for
	$(FC) $(FFLAGS) -c regridf2.for
seval.o: seval.f
	$(FC) $(FFLAGS) -c seval.f
sgf_ax_fs.o: sgf_ax_fs.f
	$(FC) $(FFLAGS) -c sgf_ax_fs.f
sgf_ax_1p.o: sgf_ax_1p.f
	$(FC) $(FFLAGS) -c sgf_ax_1p.f
#sgf_ax_str.o: sgf_ax_str.for
#	$(FC) $(FFLAGS) -c sgf_ax_str.for
#sgf_ax_s1p.o: sgf_ax_s1p.for
#	$(FC) $(FFLAGS) -c sgf_ax_s1p.for
sgf_ax_1p_ct.o: sgf_ax_1p_ct.f
	$(FC) $(FFLAGS) -c sgf_ax_1p_ct.f
#sgf_ax_gpc_ct.o: sgf_ax_gpc_ct.for
#	$(FC) $(FFLAGS) -c sgf_ax_gpc_ct.for
splc_clm.o: splc_clm.f
	$(FC) $(FFLAGS) -c splc_clm.f
strint.o: strint.for
	$(FC) $(FFLAGS) -c strint.for
thomas.o: thomas.f
	$(FC) $(FFLAGS) -c thomas.f
#
#clean:
#\rm -rf *.o pinch
#cleaner: clean
#\rm -rf output/*
#
