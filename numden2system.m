function TF_specs = numden2system(nume, denu)
    % fonction pour évaluer les paramètres temporels et fréquentiels d'une fonction de transfer
    TF_specs.num = nume;
	TF_specs.den = denu;
	
	TF_specs.TF = tf(nume, denu);
	TF_specs.BF = feedback(TF_specs.TF, 1);
	TF_specs.racines = roots(denu);
	
    step_response = stepinfo(TF_specs.TF);
	[Gm_, Pm_, Wcg_, Wcp_] = margin(TF_specs.TF);
	[mag_, pha_, wout_] = bode(TF_specs.TF);
	
	TF_specs.tr = step_response.RiseTime;
	TF_specs.ts = step_response.SettlingTime;
	TF_specs.Mp = step_response.Overshoot;
	TF_specs.Pk = step_response.Peak;
	TF_specs.Pt = step_response.PeakTime;
	
	TF_specs.Gm  = Gm_;
	TF_specs.Pm  = Pm_;
	TF_specs.Wcg = Wcg_;
	TF_specs.Wcp = Wcp_;
	
	TF_specs.mag  = mag_;
	TF_specs.pha  = pha_;
	TF_specs.wout = wout_;
end