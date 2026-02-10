%  GENERATE_REFERENCE - Generate reference data for MNPBEM regression tests.
%    Runs the same 3 demo computations as run_regression_tests.m and saves
%    the results as .mat files in tests/reference/.

disp('=== Generating MNPBEM Reference Data ===');
disp(' ');

outdir = 'tests/reference';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

%% 1. demospecstat1
disp('Generating demospecstat1 reference...');
[sca, ext, enei] = compute_demospecstat1();
save(fullfile(outdir, 'demospecstat1_ref.mat'), 'sca', 'ext', 'enei');
disp('  Done.');

%% 2. demodipstat1
disp('Generating demodipstat1 reference...');
[tot, rad, z] = compute_demodipstat1();
save(fullfile(outdir, 'demodipstat1_ref.mat'), 'tot', 'rad', 'z');
disp('  Done.');

%% 3. demoeelsstat1
disp('Generating demoeelsstat1 reference...');
[psurf, ene] = compute_demoeelsstat1();
save(fullfile(outdir, 'demoeelsstat1_ref.mat'), 'psurf', 'ene');
disp('  Done.');

disp(' ');
disp('=== Reference data generation complete ===');


%% ---- Local functions ----

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
