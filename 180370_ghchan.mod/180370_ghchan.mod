TITLE gh channel channel
: Hodgkin - Huxley h channel


NEURON {
	SUFFIX gh
	:USEION k READ ek WRITE ik
	:USEION na READ ena WRITE ina
	:RANGE ghbar, ik, ina,htau, half, slp
	NONSPECIFIC_CURRENT i
	GLOBAL inf,eh
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	ghbar =.001 (mho/cm2) <0,1e9>
	htau = 50 (ms)
	half=-80 (mV)
	slp=8 (mV)
	:ek = -77 (mV)
	:ena = 50 (mV)
}
STATE {
	n
}
ASSIGNED {
	i (mA/cm2)
	eh (mV)
	:ik (mA/cm2)
	:ina (mA/cm2)
	inf
}

INITIAL {
	rate(v)
	n = inf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	:ik = 0.7*ghbar*n*(v - ek)
	:ina = 0.3*ghbar*n*(v - ena)
	i = ghbar*n*(v - eh)
}

DERIVATIVE states {	
	rate(v)
	n' = (inf - n)/htau
}
UNITSOFF

PROCEDURE rate(v(mV)) {	
	TABLE inf DEPEND half,slp FROM -100 TO 100 WITH 200
		inf = 1/(1+exp((v-half)/slp))
}
UNITSON
