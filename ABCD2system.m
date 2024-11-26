function FT_specs = numden2tf(A, B, C, D)
    % fonction pour évaluer les paramètres temporels et fréquentiels d'une fonction de transfer
    
	FT_specs.FT = ss(A, B, C, D);
	FT_specs.BF = feedback(FT, 1);
	FT_specs.racines = roots(den);
	
    step_response = stepinfo(FT);
	[Gm_, Pm_, Wcg_, Wcp_] = margin(FT);
	[mag_, pha_, wout_] = bode(FT);
	
	FT_specs.tr = step_response.RiseTime;
	FT_specs.ts = step_response.SettlingTime;
	FT_specs.Mp = step_response.Overshoot;
	FT_specs.Pk = step_response.Peak;
	FT_specs.Pt = step_response.PeakTime;
	
	FT_specs.Gm  = Gm_;
	FT_specs.Pm  = Pm_;
	FT_specs.Wcg = Wcg_;
	FT_specs.Wcp = Wcp_;
	
	FT_specs.mag  = mag_;
	FT_specs.pha  = pha_;
	FT_specs.wout = wout_;
end