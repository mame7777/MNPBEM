function obj = clear( obj )
%  CLEAR - Clear Green functions and auxiliary matrices.
%
%  Usage for obj = bemret :
%    obj = clear( obj )

[ obj.G1i, obj.G2i, obj.L1, obj.L2, obj.Sigma1,  ...
  obj.Delta_L, obj.Delta_U, obj.Delta_p,  ...
  obj.Sigma_L, obj.Sigma_U, obj.Sigma_p ] = deal( [] );
