package com.errorgon.sensitivitytest;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class NeyerTest {

    @Test
    void round() {

        int precision = 0;
        Neyer neyer = new Neyer("in.", 0.6, 1.4, 0.1, precision);
        double results = neyer.round(10.123456789);
        Assertions.assertEquals(10, results);

        precision = 1;
        neyer = new Neyer("in.", 0.6, 1.4, 0.1, precision);
        results = neyer.round(10.123456789);
        Assertions.assertEquals(10.1, results);

        precision = 3;
        neyer = new Neyer("in.", 0.6, 1.4, 0.1, precision);
        results = neyer.round(10.123456789);
        Assertions.assertEquals(10.123, results);

    }
}