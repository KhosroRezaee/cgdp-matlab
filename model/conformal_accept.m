function ok = conformal_accept(p_act, gate, p_hemo, cfg, desc)
%CONFORMAL_ACCEPT Final accept/reject gate (paper-style)
ok = (p_act >= gate.t_act);

% Optional hemolysis constraint
if cfg.model.use_hemo
    ok = ok && (p_hemo <= cfg.model.max_p_hemo);
end

% Physchem windows
ok = ok && (desc.charge >= cfg.windows.charge_min) && (desc.charge <= cfg.windows.charge_max);
ok = ok && (desc.gravy  >= cfg.windows.gravy_min)  && (desc.gravy  <= cfg.windows.gravy_max);
ok = ok && (desc.boman  <= cfg.windows.boman_max);
ok = ok && (desc.length >= cfg.data.min_len) && (desc.length <= cfg.data.max_len);

end
