function m_out = col_prod(m_in, v)
% m_out = col_prod(m_in, v)

m_out = m_in .* (ones(size(m_in,1),1) * ToRow(v) );
