//! MolMason protocol crate

pub fn hello() -> &'static str {
    "MolMason protocol"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
