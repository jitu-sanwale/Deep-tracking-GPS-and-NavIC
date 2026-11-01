function tropo_delay = trop_UNB3(doy,lat,alt,el)
    % a diverged navigation state (NaN / off-Earth position) must not abort the
    % run: no correction then (crashed the classical receiver in the Split-S)
    if ~all(isfinite([lat alt el])) || abs(lat) > 90
        tropo_delay = 0;  return;
    end

    [K_dry, K_wet] = Trop_Saastamoinen_UNB3_Components(doy,lat,alt);
     m_dry = Trop_Black_Eisner_Map(cosd(el)); 
	 m_wet = m_dry;
     tropo_delay  = K_dry * m_dry + K_wet * m_wet;

end