function pVec = get_Permutation(pos, qmax, nStep, qmin)

        nDim = size(nStep,2);
        pVec = zeros(nDim,1);
            
        for j = 1: nDim
            seedStep = prod(nStep(j+1:end));
            if (j == nDim)
                seedStep = 1;
            end
            i = mod(floor((pos-0.9)/seedStep), nStep(j))+1;
            pVec(j) = (i-1)*(qmax(j)-qmin(j))/(nStep(j)-1)...
                       + qmin(j);
        end
end

% subroutine get_Permutation(pos, qmax, nStep, pVec, qmin, snapExtremes)
% 
%         implicit none
% 
%         !INPUT
%         integer                        , intent(in)           :: pos;
%         double precision, dimension(1:), intent(in)           :: qmax;
%         double precision, dimension(1:), intent(in), optional :: qmin;
%         integer,          dimension(1:), intent(in)           :: nStep;
%         logical, optional, intent(in) :: snapExtremes
%         !OUTPUT
%         double precision, dimension(1:), intent(out) :: pVec;
%         !LOCAL VARIABLES
%         integer :: i, j;
%         integer :: seedStep, nDim;
%         double precision :: contrib
% 
%         nDim = size(nStep);
%         contrib = 0.0d0
% 
%         if (present(snapExtremes)) then
%             if(snapExtremes) then
%                 contrib = 1.0d0
%             end if
%         end if
% 
%         if (present(qmin)) then
%             do j = 1, nDim
%                 seedStep = product(nStep(j+1:));
%                 if (j == nDim) seedStep = 1;
%                 i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
%                 pVec(j) = (dble(i)-0.5d0-contrib/2.0d0)         &
%                           *(qmax(j)-qmin(j))/(nStep(j)-contrib) &
%                           + qmin(j);
%             end do
%         else
%             do j = 1, nDim
%                 seedStep = product(nStep(j+1:));
%                 if (j == nDim) seedStep = 1;
%                 i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
%                 pVec(j) = (dble(i)-0.5d0-contrib/2.0d0) &
%                           *(qmax(j))/(nStep(j)-contrib); !qmin = 0
%             end do
%         end if
% 
%     end subroutine get_Permutation