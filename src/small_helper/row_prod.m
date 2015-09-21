function m_out = row_prod(m_in, v)
% m_out = row_prod(m_in, v)

m_out = m_in .* (ToColumn(v) * ones(1,size(m_in,2)));
