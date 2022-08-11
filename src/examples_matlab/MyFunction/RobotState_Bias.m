classdef RobotState_Bias < handle
    %UNTITLED 此处提供此类的摘要
    %   此处提供详细说明

    properties
        R_member
        v_member
        p_member
        d_all_member %3*k  k legs left,right
        bg_member
        ba_member
        
        At_member
        X_member
        Q_member
        P_member
        
        % covriance
        encoder_cov_member
        angular_vel_cov_member
        acc_cov_member
        contact_cov_member
        bg_cov_member
        ba_cov_member

        % gravity
        g_member

        %flag
        bias_flag_member;
        legcontact_flag_member %1*k   k legs
        legcontact_last_flag_member %1*k   k legs

    end

    methods
        function obj = RobotState_Bias(R,v,p,d_all,P,bias_flag)    
            % P 15*15
            obj.R_member = R;
            obj.v_member = v;
            obj.p_member = p;
            obj.d_all_member = d_all;
            obj.d_all_member = zeros(3,2);
            obj.X_member = [R,v,p;zeros(2,3),eye(2)];
            obj.P_member = P;
            obj.bias_flag_member = bias_flag;
        
            obj.bg_member = zeros(3,1);
            obj.ba_member = zeros(3,1);


            obj.encoder_cov_member = 0.01*eye(14);
            obj.angular_vel_cov_member = 0.01*eye(3);
            obj.acc_cov_member = 0.01*eye(3);
            obj.contact_cov_member = zeros(3,3,2);
            obj.contact_cov_member(:,:,1) = 0.01*eye(3);
            obj.contact_cov_member(:,:,2) = 0.01*eye(3);
            obj.g_member=[0,0,-9.81]';
            obj.legcontact_flag_member=zeros(1,2);  %1*k   k legs
            obj.legcontact_last_flag_member=zeros(1,2);
            obj.bg_cov_member = 0.01*eye(3);
            obj.ba_cov_member = 0.01*eye(3);

            obj.At_member = [zeros(3,9),-obj.R_member,zeros(3);
                             axis2skew(obj.g_member),zeros(3,6),-axis2skew(obj.v_member)*obj.R_member, - obj.R_member;
                             zeros(3),eye(3),zeros(3),-axis2skew(obj.p_member)*obj.R_member,zeros(3);
                             zeros(6,15)]; %15*15
            obj.Q_member = blkdiag(obj.angular_vel_cov_member,obj.acc_cov_member,zeros(3),obj.bg_cov_member,obj.ba_cov_member);
        end

        function [obj] = prediction(obj,w_meas,a_meas,joint_meas,Delta_T)
            %prediction X_k=f(X_{k-1})
            w_meas=w_meas(:);
            a_meas=a_meas(:);
            
            % descreted prediction
            % bias 
            % obj.bg_member = obj.bg_member;
            % obj.ba_member = obj.ba_member;

            % Base Pose Dynamics
            obj.R_member = obj.R_member * expm(axis2skew(w_meas-obj.bg_member)*Delta_T);
            obj.v_member = obj.v_member + (obj.R_member*(a_meas-obj.ba_member) + obj.g_member)*Delta_T;
            obj.p_member = obj.p_member + obj.v_member*Delta_T + ...
                            0.5*(obj.R_member*(a_meas-obj.ba_member) + obj.g_member)*Delta_T^2;
            
            % Foot Position Dynamics
            d_off= zeros(3,2);
            for i = 1:2
                d_off(:,i) = obj.p_member + obj.R_member * obj.Kin_posi(joint_meas,i);
                obj.d_all_member(:,i) = obj.d_all_member(:,i)*obj.legcontact_flag_member(i)+...
                    (1-obj.legcontact_flag_member(i))*d_off(:,i);
            end
            %update Q and At
            obj.UpdateAtandQ(joint_meas);

            % Update P
            F=eye(size(obj.At_member))+obj.At_member*Delta_T;
            if obj.bias_flag_member
                F_Q=Delta_T*blkdiag(obj.Adjoint(obj.X_member),eye(6)) ;
            else
                F_Q=Delta_T*obj.Adjoint(obj.X_member) ;
            end
            obj.P_member=F*obj.P_member*F'+F_Q*obj.Q_member*F_Q' ;

            obj.GroupX();
        end

        function obj = Correction(obj,joint_meas)
            N_now = sum(obj.legcontact_flag_member);
            R_prime = zeros(3*N_now);
            if obj.bias_flag_member
                H = zeros(3*N_now,9+3*N_now+6);
            else
                H = zeros(3*N_now,9+3*N_now);
            end

            for i_IEKF = 1:1
            % iterative 
                k=1;
                for i = 1:2
                    if obj.legcontact_flag_member(i)
                        H(k*3-2:k*3,7:9) = -eye(3);
                        H(k*3-2:k*3,9+k*3-2:9+k*3) = eye(3);
                        R_prime(k*3-2:k*3,k*3-2:k*3) = obj.R_member * obj.Kin_Jaco(joint_meas,i)*obj.encoder_cov_member...
                            *obj.Kin_Jaco(joint_meas,i)'* obj.R_member';
                        k=k+1;
                    end
                end

                K=obj.P_member*H'/(H*obj.P_member*H'+R_prime);
                inovation = zeros(N_now*3,1);
                k=1;
                for i = 1:2
                    if obj.legcontact_flag_member(i)
                        Y = zeros(5+N_now ,1);
                        Y(1:3,1) = obj.Kin_posi(joint_meas,i);
                        Y(5) =1;
                        Y(5+k) =-1;
                        inovation(k*3-2:k*3,1) = [eye(3),zeros(3,2+N_now)] * obj.X_member *Y;
                        k=k+1;
                    end
                end

                if obj.bias_flag_member
                    delta=K*inovation;
%                     K_X=K(1:9+3*N_now,:);
%                     K_bias=K(end-5:end,:);
                    
                    obj.X_member=expm(obj.Rn2liealgebra(delta(1:end-6)))*obj.X_member;
%                     bg_ba_all_inovation=K_bias * inovation;
                    obj.bg_member = obj.bg_member + delta(end-5:end-3);
                    obj.ba_member = obj.ba_member + delta(end-2:end);
%                     obj.bg_member(3)=0;
                else
                    obj.X_member=expm(obj.Rn2liealgebra(K*inovation))*obj.X_member;
                end

                obj.SeprateX();
            end
            obj.P_member=(eye(size(obj.P_member))-K*H)*obj.P_member;
%             obj.P_member=(eye(size(obj.P_member))-K*H)*obj.P_member*(eye(size(obj.P_member))-K*H)'+...
%                 K*R_prime*K';
            % TEST for consider bg as a constant
%             obj.bg_member = 
            
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
        
        function obj = ContactUpdata(obj,joint_meas)
            % Update Varible from legcontactflag
            N_now = sum(obj.legcontact_flag_member);
            N_last = sum(obj.legcontact_last_flag_member);
            
            % Group X 
            obj.X_member = eye(5+N_now);
            obj.X_member(1:3,1:3)=obj.R_member;
            obj.X_member(1:3,4)=obj.v_member;
            obj.X_member(1:3,5)=obj.p_member;
            k=1;
            for i=1:2
                if obj.legcontact_flag_member(i)
                    if obj.legcontact_last_flag_member(i)
                        obj.X_member(1:3,5+k)=obj.d_all_member(:,i);
                    else
                        obj.X_member(1:3,5+k) = obj.p_member + obj.R_member * obj.Kin_posi(joint_meas,i);
                        obj.d_all_member(:,i) = obj.X_member(1:3,5+k);
                    end
                    k=k+1;
                end
            end

            % update Cov
            M = zeros(9+3*N_now,9+3*N_last);
            G = zeros(9+3*N_now,14);
            M(1:9,1:9)=eye(9);
            k=1;
            for i=1:2
                if obj.legcontact_flag_member(i)
                    if obj.legcontact_last_flag_member(i)
                        index_I=sum(obj.legcontact_last_flag_member(1:i));
                        M((9+k*3-2):(9+k*3),(9+index_I*3-2):(9+index_I*3)) = eye(3);
                    else
                        M((9+k*3-2):(9+k*3),7:9) = eye(3);
                        G((9+k*3-2):(9+k*3),:) = obj.R_member * obj.Kin_Jaco(joint_meas,i);
                    end
                    k=k+1;
                end
            end
            if obj.bias_flag_member
                M = blkdiag(M,eye(6));G = [G;zeros(6,14)];
            end

            obj.P_member = M*obj.P_member*M'+G*obj.encoder_cov_member*G';

            obj.legcontact_last_flag_member = obj.legcontact_flag_member;
        end
        
         function [obj] = GroupX(obj)
            % Group X 
            N_now = sum(obj.legcontact_flag_member);
            obj.X_member = eye(5+N_now);
            obj.X_member(1:3,1:3)=obj.R_member;
            obj.X_member(1:3,4)=obj.v_member;
            obj.X_member(1:3,5)=obj.p_member;
            
            k=1;
            for i=1:2
                if obj.legcontact_flag_member(i)
                    obj.X_member(1:3,5+k) = obj.d_all_member(:,i);
                    k=k+1;
                end
            end
            
        end

        function obj = SeprateX(obj)
            obj.R_member = obj.X_member(1:3,1:3);
            obj.v_member = obj.X_member(1:3,4);
            obj.p_member = obj.X_member(1:3,5);
            k=1;
            for i=1:2
                if obj.legcontact_flag_member(i)
                    obj.d_all_member(:,i)=obj.X_member(1:3,5+k);
                    k=k+1;
                end
            end
        end
        
        function [obj] = UpdateAtandQ(obj,joint_meas)
            N_now = sum(obj.legcontact_flag_member);
      
            if obj.bias_flag_member
                obj.Q_member = zeros(9+3*N_now+6);
                obj.Q_member(end-5:end-3,end-5:end-3) = obj.bg_cov_member;
                obj.Q_member(end-2:end,end-2:end) = obj.ba_cov_member;
            else
                obj.Q_member = zeros(9+3*N_now);
            end
            obj.Q_member(1:3,1:3)=obj.angular_vel_cov_member;
            obj.Q_member(4:6,4:6)=obj.acc_cov_member;
            
            if obj.bias_flag_member
                obj.At_member = zeros(9+3*N_now+6);
                obj.At_member(1:9,end-5:end)=[-obj.R_member,zeros(3);
                                          -axis2skew(obj.v_member)*obj.R_member, - obj.R_member;
                                          -axis2skew(obj.p_member)*obj.R_member,zeros(3)];
            else
                obj.At_member = zeros(9+3*N_now);
            end
            obj.At_member(1:9,1:9)=[zeros(3,9);axis2skew(obj.g_member),zeros(3,6);zeros(3),eye(3),zeros(3)];
            
            k=1;
            for i=1:2
                if obj.legcontact_flag_member(i)
                    obj.Q_member(9+k*3-2:9+k*3,9+k*3-2:9+k*3) = ...
                        obj.Kin_orie(joint_meas,i)*obj.contact_cov_member(:,:,i)*obj.Kin_orie(joint_meas,i)';
                    if obj.bias_flag_member
                        obj.At_member((9+k*3-2):(9+k*3),end-5:end-3)= -axis2skew(obj.d_all_member(:,i))*obj.R_member;
                    end
                    k=k+1;
                end
            end
        end

        function J = Kin_Jaco(~,joint_meas,index_leg)
            % index_leg ---- legcontact_flag_member, 1 2 3 4 ...
            switch index_leg
                case 1
                    J = J_VectorNav_to_LeftToeBottom(joint_meas);
                case 2
                    J = J_VectorNav_to_RightToeBottom(joint_meas);
            end
        end

        function hp = Kin_posi(~,joint_meas,index_leg)
            switch index_leg
                case 1
                    hp = p_VectorNav_to_LeftToeBottom(joint_meas);
                case 2
                    hp = p_VectorNav_to_RightToeBottom(joint_meas);
            end
        end
        
        function hR = Kin_orie(~,joint_meas,index_leg)
            switch index_leg
                case 1
                    hR = R_VectorNav_to_LeftToeBottom(joint_meas);
                case 2
                    hR = R_VectorNav_to_RightToeBottom(joint_meas);
            end
        end
    end
end