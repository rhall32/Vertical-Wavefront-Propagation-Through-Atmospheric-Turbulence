%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PhaseResidues.m calculates the phase residues for a given wrapped phase
% image. Note that by convention the positions of the phase residues are 
% marked on the top left corner of the 2 by 2 regions.
%
%   active---res4---right
%      |              |
%     res1           res3
%      |              |
%   below---res2---belowright
% Phase residues with integer multiples of 2*pi are not accounted for, but 
% these rarely occur.
% Created by B.S. Spottiswoode on 07/10/2004
% Last modified on 08/10/2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residue_charge=PhaseResidues_r1(IM_phase, IM_mask);

[rows, cols]=size(IM_phase);

%The code below is simply a vectorised representation of operations on 2 by 2
%blocks in the matrix
IM_active = IM_phase;
IM_below =zeros(rows,cols);
IM_below(1:rows-1,:) = IM_phase(2:rows,:);  % Shift up 1 pixel
IM_right =zeros(rows,cols);
IM_right(:,1:cols-1) = IM_phase(:,2:cols);  % Shift left 1 pixel
IM_belowright =zeros(rows,cols);
IM_belowright(1:rows-1,1:cols-1) = IM_phase(2:rows,2:cols);  % Shift up 1 pixel & left 1 pixel

% Wrap the phase differences as we loop around the 2 by 2 blocks
res1=mod(IM_active - IM_below + pi, 2*pi) - pi;      % Up - Down (bottom row is unreliable)
res2=mod(IM_below - IM_belowright + pi, 2*pi) - pi;  % Left- Right, shifted up (bottom row & far-right column are unreliable)
res3=mod(IM_belowright - IM_right + pi, 2*pi) - pi;  % Down - Up, shifted left (bottom row & far-right column are unreliable)
res4=mod(IM_right - IM_active + pi, 2*pi) - pi;      % Right - Left, (far-right column is unreliable)

temp_residues=res1+res2+res3+res4;              %Sum the phase differences. Positive residues appear as 2*pi, negative as -2*pi.
% Far-right column and bottom row are unreliable in temp_residues
residues=(temp_residues>=6);                    %Assign 1 to positive residue (which should equal 2*pi)
residues=residues - (temp_residues<=-6);        %Assign -1 to negative residues (which should equal -2*pi)
residues(:,cols)=0; residues(rows,:)=0;         %Zero pad the border residues
%residues(:,1)=0; residues(1,:)=0;       % Not needed.  Left column and top row are accurate
residue_charge=residues;

%residue_sum=sum(sum(abs(residues)))  % DEBUG