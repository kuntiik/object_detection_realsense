_SCALE 1000 prevod na m.

    Eigen::MatrixXf plane(3, pocet bodu);
    int ii=0;
    v_mini = 0;
    for (int v = 0; v < mini.rows; v++) {
      for (int u = 0; u < mini.cols; u++) {
        unsigned char ind_p = ind_d.at<unsigned char>(v,u);
        if (ind_p>=max_i-1 && ind_p<=max_i+1) {
          float Z = (static_cast<float>(mini_ptr[v_mini + u]))/_SCALE;
          plane.col(ii++) << ((u - r_cx) * Z / r_fx), ((v - r_cy) * Z / r_fy), Z;
        }
      }
      v_mini+=mini.cols;
    }
    Eigen::Matrix<float, 3, 1> mean = plane.rowwise().mean();
    const Eigen::Matrix3Xf points_centered = plane.colwise() - mean;
    int setting = Eigen::ComputeFullU | Eigen::ComputeThinV;
    Eigen::JacobiSVD<Eigen::Matrix3Xf> svd = points_centered.jacobiSvd(setting);
    Eigen::Vector3f normal = svd.matrixU().col(2);
    Eigen::Vector3f plane_x = svd.matrixU().col(0);
    Eigen::Vector3f plane_y = svd.matrixU().col(1);
    normal.normalize();
    plane_x.normalize();
    plane_y.normalize();
    plane_d = normal(0)*mean(0)+normal(1)*mean(1)+normal(2)*mean(2);
    if (plane_d<0) {
      normal(0)=-normal(0);
      normal(1)=-normal(1);
      normal(2)=-normal(2);
      plane_d = normal(0)*mean(0)+normal(1)*mean(1)+normal(2)*mean(2);
    }
    #ifdef _DEBUG
    cout << "Normal "<< normal<<" control ii:"<<ii<<" sum hist:"<<(hist_d[max_i-1]+hist_d[max_i]+hist_d[max_i+1])<<endl;
    #endif
    cout << "Max i "<<((int)max_i)<<", "<<(max_i*0.05)<<"m aprox:"<<plane_d<<"m altitude:"<<distance<<"m poc: "<<hist_d[max_i]<<endl;

    double alf, bet;
    Eigen::Vector3f orig;
    orig[0]=pl_a; orig[1] = pl_b; orig[2]=pl_c;
    diff_vec(orig, normal, alf, bet);
    cout << "Angle diff "<<alf<<","<<bet<<endl;

    v_mini = 0;
    for (int v = 0; v < mini.rows; v++) {
      for (int u = 0; u < mini.cols; u++) {
        float val = ((float)(mini_ptr[v_mini + u])/_SCALE);  // the same scale as original depth image
        float d = (normal(0)*(u- r_cx)/r_fx + normal(1)*(v* - r_cy) / r_fy + normal(2)) * val;
        pl_d.at<float>(v,u)=d;
      }
      v_mini+=mini.cols;
    }
