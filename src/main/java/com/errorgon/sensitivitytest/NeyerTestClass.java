//package com.errorgon.sensitivitytest;
//
//import org.junit.jupiter.api.Assertions;
//import org.junit.jupiter.api.Test;
//
//import java.util.ArrayList;
//
//class NeyerTestClass {
//
//    @Test
//    public void neyerTest2() {
//        Neyer neyer = new Neyer("inches", .6, 1.4, .1, 5);
//        Boolean resultsList [] = new Boolean[]{false, true, false, true, false, true, false, true, false, true};
//        double valueList[] = new double[]{1, 1.2, 1.1, 1.2814, 1.15, 1.28059, 1.09669, 1.23314, 1.175, 1.23689};
//
//        Neyer.Run run;
//
//        for (int i = 0; i < resultsList.length; i++) {
//            run = neyer.getRun();
//            run.setResult(resultsList[i]);
//            Assertions.assertEquals(valueList[i], run.getValue());
//            neyer.setRun(run);
//            System.out.println(run);;
//        }
//    }
//
//    @Test
//    public void neyerTest() {
//        Neyer neyer = new Neyer("inches", .6, 1.4, .1, 2);
//        Boolean resultsList [] = new Boolean[]{false, false, false, false, false, true, false, false, false, false, false, false, true, false, true, false, true, true, true, true};
//        double[] muList = new double[]    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.15, 4.28, 4.52};
//        double[] sigmaList = new double[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.10, 0.19, 0.79};
//        double valueList[] = new double[]{1, 1.2, 1.4, 1.8, 2.6, 4.2, 3.4, 3.8, 4.0, 4.1, 4.28, 4.52, 5.55, 5.24, 6.37, 6.08, 7.38, 7.09, 6.89, 6.74};
//
//        Neyer.Run run;
//
//        for (int i = 0; i < resultsList.length; i++) {
//            run = neyer.getRun();
//            run.setResult(resultsList[i]);
//            Assertions.assertEquals(valueList[i], run.getValue());
//            Assertions.assertEquals(muList[i], run.getMu());
//            Assertions.assertEquals(sigmaList[i], run.getSig());
//            neyer.setRun(run);
//            System.out.println(run);;
//        }
//    }
//
//}
