{

    vector<Vec4i> hierarchy;
    vector<vector<Point>> contours;
    vector<Point> tmp;
    vector<vector<Point>> draw_cont;
    double z_avg = 0;
    Point2f corners[4];
    Matx33f R = getR(n, Vec3f(0, 0, 1));
    morphologyEx(img.clone(), img, MORPH_OPEN, Mat::ones(3,3, CV_8U));
    findContours(img.clone(), contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_SIMPLE);
    for(int i = 0; i < contours.size(); i++){
        convexHull(contours[i], tmp);
        contours[i] = tmp;
    }
    Point3f f_out;
    //for(int i = 0; i < contours.size(); i++){
    Mat res;
    for(int i = 0; i < 1; i++){
    Mat pts(3, contours[i].size(), CV_32FC1);
        for(int j=0; j < contours[i].size() ; j++){
            f_out = p2i_to_p3f(contours[i][j], depth.at<uint16_t>(contours[i][j].y,contours[i][j].x)) ;
            pts.at<float>(0,j) = f_out.x;
            pts.at<float>(1,j) = f_out.y;
            pts.at<float>(2,j) = f_out.z;
        }
        //cout << pts.col(0).size() << "  "  << endl;
        //cout << pts << endl;
        res = Mat(R)*pts;
        vector<Point2i> pts_rot;
        //Mat lay_flat(700, 1000, CV_8UC1);
        Mat lay_flat(700, 1000, CV_8UC1, Scalar(0));
        for(int k = 0; k < contours[i].size(); k++){
            pts_rot.push_back(p3f_to_p2i((Point3f)res.col(k)));
            pts_rot[k].x += 400;
            cout << pts_rot << endl;
            lay_flat.at<uint8_t>(pts_rot[k].y, pts_rot[k].x) = 255;
            z_avg += res.at<float>(2,k);
        }

        z_avg /= contours[i].size();
        cout << "z_avg" << z_avg << endl;
        draw_cont.push_back(pts_rot);
        minAreaRect(pts_rot).points(corners);
        Point2i docasny;
        Point3f corners_fin[4];
        for(int l = 0; l < 4; l++){

            cout << corners[l] << endl;
            docasny.x = corners[l].x - 400;
            docasny.y = corners[l].y;
            //docasny.z = z_avg;
            corners_fin[l] = p2i_to_p3f(docasny, z_avg); 
            cout << corners_fin[l] << endl;
        }
        Matx34f corners_as_vec;
        Vec3f dalsi_pom;
        Point2i fakt_konec[4];
        for(int l = 0; l < 4; l++){
            dalsi_pom = R.t()*(Vec3f)corners_fin[l];
            fakt_konec[l] = p3f_to_p2i((Point3f)(dalsi_pom));
            cout << fakt_konec[l] << endl;
        }
        

        for(int v = 0; v < 4; v++){

            line(img, fakt_konec[v], fakt_konec[(v+1)%4], 255);
        }
        cout << "corners" << corners << endl;
        drawContours(lay_flat, draw_cont, -1, 255);
        imshow("layed flat", lay_flat);
        imshow("fin_corner", img);
        waitKey(0);
    cout << res << endl;
    }

    
    drawContours(img, contours, 0, 255);
    RotatedRect bb;
    Point2f vert[4];
    for(int j = 0; j < contours.size(); j++){
    for(int i = 0; i < 4;i++){
        bb = minAreaRect(contours[j]);
        bb.points(vert);
        line(orig, vert[i], vert[(i+1)%4], 255);
    }
    }
    imshow("boxes", orig);
    imshow("contours", img);
    waitKey(0);
}
