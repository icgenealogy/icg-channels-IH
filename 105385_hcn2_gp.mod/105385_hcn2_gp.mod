NEURON {
	SUFFIX hcn2_gp
	NONSPECIFIC_CURRENT i
	RANGE i, g, gbar
	GLOBAL a0, b0, ah, bh, ac, bc, aa0, ba0
	GLOBAL aa0, ba0, aah, bah, aac, bac
	GLOBAL kon, koff, b, bf, ai, gca, shift
        GLOBAL eh
}

UNITS {
	(mV)	= (millivolt)
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA)	= (milliamp)
	(S)	= (siemens)
}

PARAMETER {
	gbar    = 1		(S/cm2)
	:ehcn    = -20		(mV)
	a0      = .0009		(/ms)		: parameters for alpha and beta
	b0      = .0004		(/ms)
	ah      = -95		(mV)
	bh      = -51.7		(mV)
	ac      = -.12		(/mV)
	bc      = .12		(/mV)
	aa0     = 3e-05		(/ms)		: parameters for alphaa and betaa
	ba0     = .001		(/ms)
	aah     = -94.2		(mV)
	bah     = -35.5		(mV)
	aac     = -.075		(/mV)
	bac     = .144		(/mV)
	kon     = 30		(/mM-ms)	: cyclic AMP binding parameters
	koff    = 4.5e-05	(/ms)
	b       = 80
	bf      = 8.94
	ai	= 1e-05		(mM)		: concentration cyclic AMP
	gca     = 1				: relative conductance of the bound state
	shift   = -12		(mV)		: shift in voltage dependence
	q10v    = 4				: q10 value from Magee 1998
	q10a    = 1.5				: estimated q10 for the cAMP binding reaction
	celsius			(degC)
        eh (mV)
}

ASSIGNED {
	v	(mV)
	g	(S/cm2)
	i	(mA/cm2)
	alpha	(/ms)
	beta    (/ms)
	alphaa	(/ms)
	betaa	(/ms)
}

STATE {
	c
	cac
	o
	cao
}

INITIAL {
    SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*(o + cao*gca)
	i = g*(v-eh)
}

KINETIC kin {
	LOCAL qa
	qa = q10a^((celsius-22 (degC))/10 (degC))
	rates(v)
	~ c <-> o       (alpha, beta)
	~ c <-> cac     (kon*qa*ai/bf,koff*qa*b/bf)
	~ o <-> cao     (kon*qa*ai, koff*qa)
	~ cac <-> cao   (alphaa, betaa)
	CONSERVE c + cac + o + cao = 1
}

PROCEDURE rates(v(mV)) {
	LOCAL qv	
	qv = q10v^((celsius-22 (degC))/10 (degC))
	if (v > -200) {
		alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))
	} else {
		alpha = a0*qv / (1 + exp(-((-200)-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-((-200)-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-((-200)-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-((-200)-bah-shift)*bac))
	}
}
