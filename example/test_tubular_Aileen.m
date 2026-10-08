close all; clear all;
%% Define files
dataDir = 'G:\\Aileen\\dataonly';
filenameMask = 'Annotated_OreR_Stage14_Embryo4_23_000000.tif';
filenameImage = 'Ore_ct633_antigfp488_stage14_0005-1crop_%06d.tif';
meshDir = fullfile(dataDir, "meshes");
if ~exist(meshDir, "dir")
    mkdir(meshDir)
end

%% Read in data and make mesh
origMask=readTiff4D(fullfile(dataDir,filenameMask));
padMask=padarray(logical(origMask), [1 1 1], false, 'both');
[ny, nx, nz] = size(padMask);

dx = 0.31;
dy = 0.31;
dz = 1.00;
% dx = 1;
% dy = 1;
% dz = 1;

x = (0:nx-1) * dx;
y = (0:ny-1) * dy;
z = (0:nz-1) * dz;

[F, V] = isosurface(x, y, z, double(padMask), 0.5);
VN = isonormals(x, y, z, padMask, V);
save(fullfile(meshDir,"mesh_raw.m"),"F","V");
plyFile = fullfile(meshDir, '000000.ply');
write_ply(V, F, plyFile);
fid = fopen(plyFile, 'w');
assert(fid >= 0, 'Could not open PLY file.');

fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'element vertex %d\n', size(V,1));
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'property float nx\n');
fprintf(fid, 'property float ny\n');
fprintf(fid, 'property float nz\n');
fprintf(fid, 'element face %d\n', size(F,1));
fprintf(fid, 'property list uchar int vertex_indices\n');
fprintf(fid, 'end_header\n');

fprintf(fid, '%.9g %.9g %.9g %.9g %.9g %.9g\n', ...
    [V, VN].');
% PLY uses zero-based face indices
fprintf(fid, '3 %d %d %d\n', (F - 1).');
fclose(fid);

%% Manual selection of A, P and D

figure
patch("Faces", F, "Vertices", V, ...
    "FaceColor", [0.75 0.75 0.75], ...
    "EdgeColor", "none");

axis equal
xlabel("x")
ylabel("y")
zlabel("z")
camlight headlight
lighting gouraud
rotate3d on
view(3)

%% Replace these values
A = [26, 13, 1]; %Anterior
P = [35, 10, 23];%Posterier
D = [40, 50, 11];%Somewhere in between

%% Check
hold on
plot3(A(1), A(2), A(3), "ro", "MarkerSize", 12, "LineWidth", 2)
plot3(P(1), P(2), P(3), "bo", "MarkerSize", 12, "LineWidth", 2)
plot3(D(1), D(2), D(3), "go", "MarkerSize", 12, "LineWidth", 2)

text(A(1), A(2), A(3), "  A", "Color","r", "FontSize",14)
text(P(1), P(2), P(3), "  P", "Color","b", "FontSize",14)
text(D(1), D(2), D(3), "  D", "Color","g", "FontSize",14)

legend("mesh","A","P","D")
save("APD_points.mat","A","P","D");
% Comment line below to view figure
%close all

%% Minimal inputs for TubULAR object
xp=struct();
xp.fileMeta = struct();
xp.fileMeta.dataDir = dataDir;
xp.fileMeta.filenameFormat = filenameImage;
xp.fileMeta.timePoints = 0;
xp.fileMeta.stackResolution = [0.31, 0.31, 1.00];
xp.fileMeta.nChannels = 1;
xp.fileMeta.swapZT = false;
xp.fileMeta.timeStampStringSpec = '%06d';

xp.expMeta = struct();
xp.expMeta.channelsUsed = 1;

xp.detectOptions = struct();
xp.detectOptions.ssfactor=1;
xp.detectOptions.ofn_smoothply = '';

opts=struct();
opts.flipy = false;
opts.meshDir = meshDir;
opts.timeUnits = 'frame';
opts.spaceUnits = "um";
% Resolution of rectangular pullback, start low and increase
opts.nU = 200;
opts.nV = 200;
opts.lambda_err = 0;
opts.nmodes     = 0;
opts.zwidth     = 0;

tubi = TubULAR(xp,opts);

rawMesh.v=V;
rawMesh.f=F;
tubi.currentMesh.rawMesh=rawMesh;
tubi.setTime(0);

apdvCoordsOpts = struct();
% change to true to select
apdvCoordsOpts.overwrite = false;
apdvCoordsOpts.preview = false;

tubi.computeAPDVCoords(apdvCoordsOpts);

apdvOpts = struct();
apdvOpts.timePoints = 0;
apdvOpts.use_iLastik = false;
% change to true to select
apdvOpts.overwrite = false;
apdvOpts.preview = false;

%
tubi.computeAPDpoints(apdvOpts);

meshAlignOpts=struct();
meshAlignOpts.overwrite=true;
meshAlignOpts.normal_step=0;

tubi.alignMeshesAPDV(meshAlignOpts);
%haven't tested this yet
%tubi.plotSeriesOnSurfaceTexturePatch() ;

cntrlineOpts=struct();
cntrlineOpts.timePoints=0;
cntrlineOpts.dilation=0;

tubi.generateFastMarchingCenterlines(cntrlineOpts);

sliceOpts=struct();

methodOpts = struct();
methodOpts.overwrite = true;

tubi.sliceMeshEndcaps(sliceOpts, methodOpts);
tubi.cleanCylMeshes();

    tubi.setTime(0) ;

    % Create the Cut Mesh
    tubi.generateCurrentCutMesh([])

        % Plot the cutPath (cutP) in 3D
    tubi.plotCutPath(tubi.currentMesh.cutMesh, tubi.currentMesh.cutPath)

        tubi.getCurrentUVCutMesh() ;

    spcutMeshOptions = struct() ;
    spcutMeshOptions.t0_for_phi0 = tubi.t0set() ;  % which timepoint do we define corners of pullback map
    spcutMeshOptions.save_phi0patch = false ;
    spcutMeshOptions.iterative_phi0 = false ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    tubi.generateCurrentSPCutMesh([], spcutMeshOptions) ;

    pbOptions = struct() ;
    pbOptions.numLayers = [0 1] ; % how many onion layers over which to take MIP
    tubi.imSize = 400;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;