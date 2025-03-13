CREATE DATABASE IF NOT EXISTS gene_set_enrichment;
USE gene_set_enrichment;

-- Gene Table: gid, gname. May add [attributes]
CREATE TABLE IF NOT EXISTS Genes(
    gid INT AUTO_INCREMENT PRIMARY KEY,
    gname VARCHAR(255) NOT NULL,
    UNIQUE(gname)
);

-- Pathway Table: pid, size, weight. May add [source]
CREATE TABLE IF NOT EXISTS Pathways(
    pid INT AUTO_INCREMENT PRIMARY KEY,
    pname VARCHAR(255) NOT NULL,
    size INT,
    weight DOUBLE
);

-- Gene_Pathway Table: (pid, gid)
CREATE TABLE IF NOT EXISTS Gene_Pathway(
    pid INT,
    gid INT,
    PRIMARY KEY(pid, gid),
    FOREIGN KEY (pid) REFERENCES Pathways (pid) ON DELETE CASCADE,
    FOREIGN KEY (gid) REFERENCES Genes (gid) ON DELETE CASCADE
);

/*
Files Table: (fname, scale)
CONDITIONS:
1. fname is PRIMARY KEY.
2. scale should be an double >= 1. DEFAULT value is 1
*/
CREATE TABLE IF NOT EXISTS Files(
    fname VARCHAR(255) PRIMARY KEY,
    scale DOUBLE DEFAULT 1
);

DELIMITER $$

CREATE TRIGGER scale_insert
BEFORE INSERT ON Files
FOR EACH ROW
BEGIN
    DECLARE msg TEXT;
    IF NEW.scale < 1 THEN
        SET msg = CONCAT('Requirement violated: ', NEW.fname, ' must have an optional scale >= 1.');
        SIGNAL SQLSTATE '45000'
        SET MESSAGE_TEXT = msg;
    END IF;
END$$

DELIMITER ;

/*
File_Pathway Table: (fname, pid)
# Assume that one file will not contain multiple pathways with the same name
# However, different files may have duplicate pathway(name)s
*/
CREATE TABLE IF NOT EXISTS File_Pathway(
    fname VARCHAR(255),
    pid INT,
    PRIMARY KEY (fname, pid),
    FOREIGN KEY (fname) REFERENCES Files (fname) ON DELETE CASCADE,
    FOREIGN KEY (pid) REFERENCES Pathways (pid) ON DELETE CASCADE
);