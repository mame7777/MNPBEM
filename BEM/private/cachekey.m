function key = cachekey( enei )
%CACHEKEY Build cache key for wavelength-dependent data.
%  INPUT
%    enei : wavelength (nm) used to tag cached preconditioners (caller ensures units)
%  OUTPUT
%    key  : valid struct field name derived from ENEI
precision = 15;
%  keep enough digits to distinguish double precision wavelengths.
fmt = sprintf( '%%0.%df', precision );
key = matlab.lang.makeValidName( [ 'e_', sprintf( fmt, enei ) ] );
