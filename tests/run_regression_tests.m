%  RUN_REGRESSION_TESTS - Run regression tests for MNPBEM.
%    Computes results from 3 demo simulations and compares them against
%    reference data stored in tests/reference/*.mat.
%    Asserts relative error within 1e-6 for all quantities.

disp('=== MNPBEM Regression Tests ===');
disp(' ');

tol = 1e-6;
all_passed = true;

%% Test 1: demospecstat1 - Gold nanosphere scattering
disp('--- Test 1: demospecstat1 (gold nanosphere scattering) ---');
try
    [sca, ext, enei] = compute_demospecstat1();
    ref = load('tests/reference/demospecstat1_ref.mat');

    assert_close(sca, ref.sca, tol, 'demospecstat1: sca');
    assert_close(ext, ref.ext, tol, 'demospecstat1: ext');
    disp('PASSED');
catch e
    disp(['FAILED: ' e.message]);
    all_passed = false;
end
disp(' ');

%% Test 2: demodipstat1 - Dipole decay rate
disp('--- Test 2: demodipstat1 (dipole decay rate) ---');
try
    [tot, rad, z] = compute_demodipstat1();
    ref = load('tests/reference/demodipstat1_ref.mat');

    assert_close(tot, ref.tot, tol, 'demodipstat1: tot');
    assert_close(rad, ref.rad, tol, 'demodipstat1: rad');
    disp('PASSED');
catch e
    disp(['FAILED: ' e.message]);
    all_passed = false;
end
disp(' ');

%% Test 3: demoeelsstat1 - EELS loss probability
disp('--- Test 3: demoeelsstat1 (EELS loss probability) ---');
try
    [psurf, ene] = compute_demoeelsstat1();
    ref = load('tests/reference/demoeelsstat1_ref.mat');

    assert_close(psurf, ref.psurf, tol, 'demoeelsstat1: psurf');
    disp('PASSED');
catch e
    disp(['FAILED: ' e.message]);
    all_passed = false;
end
disp(' ');

%% Summary
if all_passed
    disp('=== All regression tests PASSED ===');
else
    error('Some regression tests FAILED.');
end


%% ---- Local functions ----

function assert_close(actual, expected, tol, name)
%  ASSERT_CLOSE - Assert relative error is within tolerance.
    denom = max(abs(expected(:)));
    if denom == 0
        rel_err = max(abs(actual(:) - expected(:)));
    else
        rel_err = max(abs(actual(:) - expected(:))) / denom;
    end
    if rel_err > tol
        error('%s: relative error %e exceeds tolerance %e', name, rel_err, tol);
    end
    fprintf('  %s: max relative error = %e (OK)\n', name, rel_err);
end

function [sca, ext, enei] = compute_demospecstat1()
%  Compute gold nanosphere scattering cross sections (from demospecstat1).
    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv' );
    epstab = { epsconst( 1 ), epstable( 'gold.dat' ) };
    diameter = 10;
    p = comparticle( epstab, { trisphere( 144, diameter ) }, [ 2, 1 ], 1, op );

    bem = bemsolver( p, op );
    exc = planewave( [ 1, 0, 0; 0, 1, 0 ], [ 0, 0, 1; 0, 0, 1 ], op );
    enei = linspace( 400, 700, 80 );
    sca = zeros( length( enei ), 2 );
    ext = zeros( length( enei ), 2 );

    for ien = 1 : length( enei )
        sig = bem \ exc( p, enei( ien ) );
        sca( ien, : ) = exc.sca( sig );
        ext( ien, : ) = exc.ext( sig );
    end
end

function [tot, rad, z] = compute_demodipstat1()
%  Compute dipole decay rates above gold nanosphere (from demodipstat1).
    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv' );
    epstab = { epsconst( 1 ), epstable( 'gold.dat' ) };
    diameter = 10;
    p = comparticle( epstab, { trisphere( 144, diameter ) }, [ 2, 1 ], 1, op );

    enei = 550;
    z = reshape( linspace( 0.6, 1.5, 51 ) * diameter, [], 1 );
    pt = compoint( p, [ 0 * z, 0 * z, z ] );
    dip = dipole( pt, [ 1, 0, 0; 0, 0, 1 ], op );

    bem = bemsolver( p, op );
    sig = bem \ dip( p, enei );
    [ tot, rad ] = dip.decayrate( sig );
end

function [psurf, ene] = compute_demoeelsstat1()
%  Compute EELS loss probability for silver nanosphere (from demoeelsstat1).
    eV2nm = 1 / 8.0655477e-4;

    op = bemoptions( 'sim', 'stat', 'waitbar', 0, 'interp', 'curv' );
    epstab = { epsconst( 1 ), epstable( 'silver.dat' ) };
    diameter = 30;
    p = comparticle( epstab, { trisphere( 144, diameter ) }, [ 2, 1 ], 1, op );

    [ width, vel ] = deal( 0.5, eelsbase.ene2vel( 200e3 ) );
    imp = 5;
    ene = linspace( 3, 4.5, 40 );
    enei = eV2nm ./ ene;

    bem = bemsolver( p, op );
    exc = electronbeam( p, [ diameter / 2 + imp, 0 ], width, vel, op );
    psurf = zeros( size( ene ) );

    for ien = 1 : length( enei )
        sig = bem \ exc( enei( ien ) );
        psurf( ien ) = exc.loss( sig );
    end
end
