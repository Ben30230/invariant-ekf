classdef RobotState < handle
    %UNTITLED 此处提供此类的摘要
    %   此处提供详细说明

    properties
        R_member
        v_member
        p_member
        dl_member
        dr_member
        X_member
        P_member

        % covriance
        encoder_cov_member
        angular_vel_cov_member
        acc_cov_member
        vel_leftcontact_cov_member
        vel_rightcontact_cov_member
        
        % gravity
        g_member

        %flag
        contact_flag_member %0:none 1: left only 2: right only 3: Two 
    end

    methods
        function obj = RobotState(R,v,p,dl,dr,P)    
            % P 15*15
            obj.R_member = R;
            obj.v_member = v;
            obj.p_member = p;
            obj.dl_member = dl;
            obj.dr_member = dr;
            obj.X_member = [R,v,p,dl,dr;zeros(4,3),eye(4)];
            obj.P_member = P;
            
            if nargin < 7
                obj.encoder_cov_member = 0.01*eye(14);
                obj.angular_vel_cov_member = 0.01*eye(3);
                obj.acc_cov_member = 0.01*eye(3);
                obj.vel_leftcontact_cov_member = 0.01*eye(3);
                obj.vel_rightcontact_cov_member = 0.01*eye(3);
                obj.g_member=[0,0,-9.81]';
                obj.contact_flag_member=3;
            end
        end

        function [obj] = prediction(obj,w_meas,a_meas,joint_meas,Delta_T)
            %prediction X_k=f(X_{k-1})
            w_meas=w_meas(:);
            a_meas=a_meas(:);
            switch obj.contact_flag_member
                case 0
                    At=[zeros(3,9);axis2skew(obj.g_member),zeros(3,6);zeros(3),eye(3),zeros(3)];  % 9*9
                    leftcontact_flag=0;
                    rightcontact_flag=0;
                case 1
                    At=[zeros(3,12);axis2skew(obj.g_member),zeros(3,9);zeros(3),eye(3),zeros(3,6)
                        zeros(3,12)];%12*12
                    leftcontact_flag=1;
                    rightcontact_flag=0;
                case 2
                    At=[zeros(3,12);axis2skew(obj.g_member),zeros(3,9);zeros(3),eye(3),zeros(3,6)
                        zeros(3,12)];%12*12
                    leftcontact_flag=0;
                    rightcontact_flag=1;
                case 3
                    At=[zeros(3,15);axis2skew(obj.g_member),zeros(3,12);zeros(3),eye(3),zeros(3,9)
                        zeros(6,15)];%15*15
                    leftcontact_flag=1;
                    rightcontact_flag=1;
            end
            % Generate X-matrix form
            obj.X_member(1:3,1:3)=obj.R_member;
            obj.X_member(1:3,4)=obj.v_member;
            obj.X_member(1:3,5)=obj.p_member;
            % cov
            h_R_left=R_VectorNav_to_LeftToeBottom(joint_meas);
            h_R_right=R_VectorNav_to_RightToeBottom(joint_meas);
            Q=zeros(3*3+3*(leftcontact_flag+rightcontact_flag));
            Q(1:3,1:3)=obj.angular_vel_cov_member;
            Q(4:6,4:6)=obj.acc_cov_member;
            if leftcontact_flag && rightcontact_flag
                obj.X_member(1:3,6)=obj.dl_member;
                obj.X_member(1:3,7)=obj.dr_member;
                Q(10:12,10:12)=h_R_left*obj.vel_leftcontact_cov_member*h_R_left';
                Q(13:15,13:15)=h_R_right*obj.vel_rightcontact_cov_member*h_R_right';
            elseif leftcontact_flag
                obj.X_member(1:3,6)=obj.dl_member;
                Q(10:12,10:12)=h_R_left*obj.vel_leftcontact_cov_member*h_R_left';
            elseif rightcontact_flag
                obj.X_member(1:3,6)=obj.dr_member;
                Q(10:12,10:12)=h_R_right*obj.vel_rightcontact_cov_member*h_R_right';
            end

            F=eye(size(At))+At*Delta_T;
            F_Q=Delta_T*obj.Adjoint(obj.X_member);

            % descreted prediction
            % Base Pose Dynamics
            obj.R_member = obj.R_member * expm(axis2skew(w_meas*Delta_T));
            obj.v_member = obj.v_member + (obj.R_member*a_meas + obj.g_member)*Delta_T;
            obj.p_member = obj.p_member + obj.v_member*Delta_T + ...
                            0.5*(obj.R_member*a_meas + obj.g_member)*Delta_T^2;
            
            % Foot Position Dynamics
            dL_off = obj.p_member + obj.R_member * p_VectorNav_to_LeftToeBottom(joint_meas);  % {W}_p_{WL}
            dR_off = obj.p_member + obj.R_member * p_VectorNav_to_RightToeBottom(joint_meas); % {W}_p_{WR}
            obj.dl_member = leftcontact_flag*obj.dl_member + (1-leftcontact_flag)*dL_off;
            obj.dr_member = rightcontact_flag*obj.dr_member + (1-rightcontact_flag)*dR_off;

            obj.P_member=F*obj.P_member*F'+F_Q*Q*F_Q' ;
            
            obj.StateGroup();
        end

        function obj = Correction(obj,joint_meas)
            switch obj.contact_flag_member
                case 0
                    
                case 1
                    H=[zeros(3,6),-eye(3),eye(3)];  % 3*12
                    Jacobian_left=J_VectorNav_to_LeftToeBottom(joint_meas);
                    R_prime=obj.R_member*Jacobian_left*obj.encoder_cov_member*Jacobian_left'*obj.R_member';
                    K=obj.P_member*H'/(H*obj.P_member*H'+R_prime);
                    inovation=[eye(3),zeros(3)]*obj.X_member*[p_VectorNav_to_LeftToeBottom(joint_meas);0;1;-1];%3*1
                    obj.X_member=expm(obj.Rn2liealgebra(K*inovation))*obj.X_member;
                    obj.P_member=(eye(12)-K*H)*obj.P_member;
                case 2
                    H=[zeros(3,6),-eye(3),eye(3)];  % 3*12
                    Jacobian_right=J_VectorNav_to_RightToeBottom(joint_meas);
                    R_prime=obj.R_member*Jacobian_right*obj.encoder_cov_member*Jacobian_right'*obj.R_member';
                    K=obj.P_member*H'/(H*obj.P_member*H'+R_prime);
                    inovation=[eye(3),zeros(3)]*obj.X_member*[p_VectorNav_to_RightToeBottom(joint_meas);0;1;-1];%3*1
                    obj.X_member=expm(obj.Rn2liealgebra(K*inovation))*obj.X_member;
                    obj.P_member=(eye(12)-K*H)*obj.P_member;
                case 3
                    H=[zeros(3,6),-eye(3),eye(3),zeros(3);
                       zeros(3,6),-eye(3),zeros(3),eye(3);];  % 6*15
                    Jacobian_left=J_VectorNav_to_LeftToeBottom(joint_meas);
                    Jacobian_right=J_VectorNav_to_RightToeBottom(joint_meas);
                    R_prime_left=obj.R_member*Jacobian_left*obj.encoder_cov_member*Jacobian_left'*obj.R_member';
                    R_prime_right=obj.R_member*Jacobian_right*obj.encoder_cov_member*Jacobian_right'*obj.R_member';
                    K=obj.P_member*H'/(H*obj.P_member*H'+blkdiag(R_prime_left,R_prime_right));%15*6
                    inovation_left=[eye(3),zeros(3,4)]*obj.X_member*[p_VectorNav_to_LeftToeBottom(joint_meas);0;1;-1;0];%3*1
                    inovation_right=[eye(3),zeros(3,4)]*obj.X_member*[p_VectorNav_to_RightToeBottom(joint_meas);0;1;0;-1];%3*1
                    inovation_all=[inovation_left;inovation_right];%6*1
                    obj.X_member=expm(obj.Rn2liealgebra(K*inovation_all))*obj.X_member;
                    obj.P_member=(eye(15)-K*H)*obj.P_member;
            end
            obj.StateSperate;
        end

        function AdjX = Adjoint(~,X)
            % Adjoint of SE_N(3)         
            N = size(X,2)-3;
            R = X(1:3,1:3);
            R_cell = repmat({R}, 1, N+1); 
            AdjX = blkdiag(R_cell{:});
            for i = 1:N
                AdjX(3*i+1:3*i+3,1:3) = axis2skew(X(1:3,i+3))*R;
            end
        end

        function output = Rn2liealgebra(~,xi)
            xi=xi(:);
            N=length(xi)/3+2;
            output=zeros(N);
            for i=3:N
                if i==3
                    output(1:3,1:3)=axis2skew(xi(1:3));
                else
                    output(1:3,i)=xi(((i-3)*3+1):((i-3)*3+3));
                end
            end
        end

        function obj = Contactnone2left(obj,joint_meas)
            % change X
            obj.contact_flag_member = 1;
            X_temp=obj.X_member;
            obj.X_member = zeros(6);
            obj.X_member(1:5,1:5)=X_temp;
            obj.dl_member=obj.p_member+obj.R_member*p_VectorNav_to_LeftToeBottom(joint_meas);
            obj.X_member(1:3,6)=obj.dl_member;
            obj.X_member(4:6,4:6)=eye(3);
            
            % change Cov
            M=[eye(9);zeros(3,6),eye(3)];   %12*9
            Jacobian_left=J_VectorNav_to_LeftToeBottom(joint_meas);
            G=[zeros(9,14);obj.R_member*Jacobian_left];
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';
        end

        function obj = Contactnone2right(obj,joint_meas)
            % change X
            obj.contact_flag_member = 2;   
            X_temp=obj.X_member;
            obj.X_member = zeros(6);
            obj.X_member(1:5,1:5)=X_temp;
            obj.dr_member=obj.p_member+obj.R_member*p_VectorNav_to_RightToeBottom(joint_meas);
            obj.X_member(1:3,6)=obj.dr_member;
            obj.X_member(4:6,4:6)=eye(3);
            
            % change Cov
            M=[eye(9);zeros(3,6),eye(3)];   %12*9
            Jacobian_right=J_VectorNav_to_RightToeBottom(joint_meas);
            G=[zeros(9,14);obj.R_member*Jacobian_right];
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';
        end

        function obj = Contactnone2two(obj,joint_meas)
            % change X
            obj.contact_flag_member = 3;
            X_temp=obj.X_member;
            obj.X_member = zeros(7);
            obj.X_member(1:5,1:5)=X_temp;
            obj.dl_member=obj.p_member+obj.R_member*p_VectorNav_to_LeftToeBottom(joint_meas);
            obj.dr_member=obj.p_member+obj.R_member*p_VectorNav_to_RightToeBottom(joint_meas);
            obj.X_member(1:3,6)=obj.dl_member;
            obj.X_member(1:3,7)=obj.dr_member;
            obj.X_member(4:7,4:7)=eye(4);
            
            % change Cov
            M=[eye(9);zeros(3,6),eye(3);zeros(3,6),eye(3)];   %15*9
            Jacobian_left=J_VectorNav_to_LeftToeBottom(joint_meas);
            Jacobian_right=J_VectorNav_to_RightToeBottom(joint_meas);
            G=[zeros(9,14);obj.R_member*Jacobian_left;obj.R_member*Jacobian_right]; %15*14
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';
        end
        
        function obj = Contactleft2right(obj,joint_meas)
            % change X
            obj.contact_flag_member = 2;
            obj.dr_member=obj.p_member+obj.R_member*p_VectorNav_to_RightToeBottom(joint_meas);
            obj.X_member(1:3,6)=obj.dr_member;
            
            % change Cov
            M=[eye(9),zeros(9,3);zeros(3,6),eye(3),zeros(3)];   %12*12
            Jacobian_right=J_VectorNav_to_RightToeBottom(joint_meas);
            G=[zeros(9,14);obj.R_member*Jacobian_right]; %12*14
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';

        end
        
        function obj = Contactleft2none(obj,~)
            % change X
            obj.contact_flag_member = 0;
            obj.X_member=obj.X_member(1:5,1:5);
            
            % change Cov
            M=[eye(9),zeros(9,3)];   %9*12
            obj.P_member = M *obj.P_member *M';
        end
        
        function obj = Contactleft2two(obj,joint_meas)
            % change X
            obj.contact_flag_member = 3;
            X_temp=obj.X_member;
            obj.X_member = zeros(7);
            obj.X_member(1:6,1:6)=X_temp;
            obj.dr_member=obj.p_member+obj.R_member*p_VectorNav_to_RightToeBottom(joint_meas);
            obj.X_member(1:3,7)=obj.dr_member;
            obj.X_member(end,end)=1;
            
            % change Cov
            M=[eye(12);zeros(3,6),eye(3),zeros(3)];   %15*12
            Jacobian_right=J_VectorNav_to_RightToeBottom(joint_meas);
            G=[zeros(12,14);obj.R_member*Jacobian_right]; %15*14
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';
        end
       
        function obj = Contactright2left(obj,joint_meas)
            % change X
            obj.contact_flag_member = 1;
            obj.dl_member=obj.p_member+obj.R_member*p_VectorNav_to_LeftToeBottom(joint_meas);
            obj.X_member(1:3,6)=obj.dl_member;
            
            % change Cov
            M=[eye(9),zeros(9,3);zeros(3,6),eye(3),zeros(3)];   %12*12
            Jacobian_left=J_VectorNav_to_LeftToeBottom(joint_meas);
            G=[zeros(9,14);obj.R_member*Jacobian_left]; %12*14
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';
        end
        
        function obj = Contactright2none(obj,~)
            % change X
            obj.contact_flag_member = 0;
            obj.X_member=obj.X_member(1:5,1:5);
            
            % change Cov
            M=[eye(9),zeros(9,3)];   %9*12
            obj.P_member = M *obj.P_member *M';
        end
      
        function obj = Contactright2two(obj,joint_meas)
            % change X
            obj.contact_flag_member = 3;
            X_temp=obj.X_member;
            obj.X_member = zeros(7);
            obj.X_member(1:6,1:6)=X_temp;
            obj.dl_member=obj.p_member+obj.R_member*p_VectorNav_to_LeftToeBottom(joint_meas);
            obj.X_member(1:3,6)=obj.dl_member;
            obj.X_member(1:3,7)=X_temp(1:3,6);
            obj.X_member(end,end)=1;
            
            % change Cov
            M=[eye(9),zeros(9,3);zeros(3,6),eye(3),zeros(3);zeros(3,9),eye(3)];   %15*12
            Jacobian_left=J_VectorNav_to_LeftToeBottom(joint_meas);
            G=[zeros(9,14);obj.R_member*Jacobian_left;zeros(3,14)]; %15*14
            obj.P_member = M *obj.P_member *M' + G*obj.encoder_cov_member*G';
        end

        function obj = Contacttwo2left(obj,~)
            % change X
            obj.contact_flag_member = 1;
            obj.X_member = obj.X_member(1:6,1:6);
            
            % change Cov
            M=[eye(12),zeros(12,3)];   %12*15
            obj.P_member = M *obj.P_member *M';
        end

        function obj = Contacttwo2right(obj,~)
            % change X
            obj.contact_flag_member = 2;
            X_temp=obj.X_member;
            obj.X_member = obj.X_member(1:6,1:6);
            obj.X_member(1:3,6)=X_temp(1:3,end);
            
            % change Cov
            M=[eye(9),zeros(9,6);zeros(3,12),eye(3)];   %12*15
            obj.P_member = M *obj.P_member *M';
        end

        function obj = Contacttwo2none(obj,~)
            % change X
            obj.contact_flag_member = 0;
            obj.X_member = obj.X_member(1:5,1:5);
            
            % change Cov
            M=[eye(9),zeros(9,6)];   %9*15
            obj.P_member = M *obj.P_member *M';
        end

        function obj = StateSperate(obj)
            obj.R_member = obj.X_member(1:3,1:3);
            obj.v_member = obj.X_member(1:3,4);
            obj.p_member = obj.X_member(1:3,5);

            if obj.contact_flag_member == 1
                obj.dl_member = obj.X_member(1:3,6);
            elseif obj.contact_flag_member == 2
                obj.dr_member = obj.X_member(1:3,6);
            elseif obj.contact_flag_member == 3
                obj.dl_member = obj.X_member(1:3,6);
                obj.dr_member = obj.X_member(1:3,7);
            end
        end
        
        function [obj] = StateGroup(obj)
            obj.X_member(1:3,1:3) = obj.R_member;
            obj.X_member(1:3,4) = obj.v_member;
            obj.X_member(1:3,5) = obj.p_member;

            if obj.contact_flag_member == 1
                obj.X_member(1:3,6) = obj.dl_member;
            elseif obj.contact_flag_member == 2
                obj.X_member(1:3,6) = obj.dr_member;
            elseif obj.contact_flag_member == 3
                obj.X_member(1:3,6) = obj.dl_member;
                obj.X_member(1:3,7) = obj.dr_member;
            end
        end

    end
end